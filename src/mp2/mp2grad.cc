//
// BAGEL - Parallel electron correlation program.
// Filename: mp2grad.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <stddef.h>
#include <src/mp2/mp2grad.h>
#include <src/grad/cphf.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <src/util/f77.h>
#include <src/smith/prim_op.h>
#include <src/prop/dipole.h>
#include <src/grad/gradeval.h>

using namespace std;
using namespace bagel;

MP2Grad::MP2Grad(const multimap<string, string> input, const shared_ptr<const Geometry> g) : MP2(input, g, shared_ptr<const Reference>()) {

}


// does nothing.
void MP2Grad::compute() { }


template<>
shared_ptr<GradFile> GradEval<MP2Grad>::compute() {
  auto tp1 = chrono::high_resolution_clock::now();

  const size_t ncore = ref_->ncore();

  // since this is only for closed shell
  const size_t naux = geom_->naux();
  const size_t nocc = ref_->nocc() - ncore;
  const size_t nocca = ref_->nocc();
  if (nocc < 1) throw runtime_error("no correlated electrons");
  const size_t nvirt = geom_->nbasis() - nocca;
  if (nvirt < 1) throw runtime_error("no virtuals orbitals");
  assert(geom_->nbasis() == ref_->coeff()->mdim());

  const size_t nbasis = geom_->nbasis();

  const double* const coeff = ref_->coeff()->data() + ncore*nbasis;
  const double* const ocoeff = ref_->coeff()->data();
  const double* const vcoeff = coeff + nocc*nbasis;

  // first compute half transformed integrals
  shared_ptr<const DF_Half> half = geom_->df()->compute_half_transform(coeff, nocc);
  // TODO this is a waste...
  shared_ptr<const DF_Half> halfjj = geom_->df()->compute_half_transform(ocoeff, nocca)->apply_JJ();
  // second transform for virtual index
  // this is now (naux, nocc, nvirt)
  shared_ptr<const DF_Full> full = half->compute_second_transform(vcoeff, nvirt)->apply_J();
  shared_ptr<const DF_Full> bv = full->apply_J();
  shared_ptr<DF_Full> gia = bv->clone();

  auto tp2 = chrono::high_resolution_clock::now();
  auto dr1 = chrono::duration_cast<chrono::milliseconds>(tp2-tp1);
  cout << setw(60) << left << "    * 3-index integral transformation done" << right << setw(10) << setprecision(2) << dr1.count()*0.001 << endl << endl;

  // assemble
  unique_ptr<double[]> buf(new double[nocc*nvirt*nocc]); // it is implicitly assumed that o^2v can be kept in core in each node
  unique_ptr<double[]> buf2(new double[nocc*nvirt*nocc]);
  vector<double> eig_tm = ref_->eig();
  vector<double> eig(eig_tm.begin()+ncore, eig_tm.end());

  shared_ptr<Matrix> dmp2(new Matrix(nbasis, nbasis));
  double* optr = dmp2->element_ptr(ncore, ncore);
  double* vptr = dmp2->element_ptr(nocca, nocca);

  double ecorr = 0.0;
  for (size_t i = 0; i != nvirt; ++i) {
    // nocc * nvirt * nocc
    unique_ptr<double[]> data = full->form_4index(full, i);
    copy(data.get(), data.get()+nocc*nvirt*nocc, buf.get());
    copy(data.get(), data.get()+nocc*nvirt*nocc, buf2.get());

    // using SMITH's symmetrizer (src/smith/prim_op.h)
    SMITH::sort_indices<2,1,0,2,1,-1,1>(data, buf, nocc, nvirt, nocc);
    double* tdata = buf.get();
    double* bdata = buf2.get();
    for (size_t j = 0; j != nocc; ++j) {
      for (size_t k = 0; k != nvirt; ++k) {
        for (size_t l = 0; l != nocc; ++l, ++tdata, ++bdata) {
          const double denom = 1.0 / (-eig[i+nocc]+eig[j]-eig[k+nocc]+eig[l]);
          *tdata *= denom;
          *bdata *= denom;
        }
      }
    }
    ecorr += ddot_(nocc*nvirt*nocc, data, 1, buf, 1);

    // form Gia : TODO distribute
    // Gia(D|ic) = BV(D|ja) G_c(ja|i)
    // BV and gia are DF_Full 
    const size_t offset = i*nocc*naux;
    gia->set_product(bv, buf, nocc, offset);

    // G(ja|ic) -> G_c(a,ij)
    SMITH::sort_indices<1,2,0,0,1,1,1>(buf, data, nocc, nvirt, nocc);
    // T(jb|ic) -> T_c(b,ij)
    SMITH::sort_indices<1,2,0,0,1,1,1>(buf2, buf, nocc, nvirt, nocc);
    // D_ab = G(ja|ic) T(jb|ic)
    dgemm_("N", "T", nvirt, nvirt, nocc*nocc, 2.0, buf.get(), nvirt, data.get(), nvirt, 1.0, vptr, nbasis);
    // D_ij = - G(ja|kc) T(ia|kc)
    dgemm_("T", "N", nocc, nocc, nvirt*nocc, -2.0, buf.get(), nvirt*nocc, data.get(), nvirt*nocc, 1.0, optr, nbasis);
  }

  auto tp3 = chrono::high_resolution_clock::now();
  auto dr2 = chrono::duration_cast<chrono::milliseconds>(tp3-tp2);
  cout << left << setw(60) << "    * assembly (+ unrelaxed density matrices) done" << right << setw(10) << setprecision(2) << dr2.count()*0.001 << endl << endl;
  cout << "      MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << ecorr << endl << endl;

  // L''aq = 2 Gia(D|ia) (D|iq)
  unique_ptr<double[]> lai(new double[nocca*nvirt]);
  const unique_ptr<double[]> laq = gia->form_2index(half, 2.0);
  {
    dgemm_("N", "N", nvirt, nocca, nbasis, 1.0, laq.get(), nvirt, ocoeff, nbasis, 0.0, lai.get(), nvirt);
  }

  // Gip = Gia(D|ia) C+_ap
  shared_ptr<DF_Half> gip = gia->back_transform(vcoeff);
  // Liq = 2 Gip(D) * (D|pq)
  unique_ptr<double[]> lia(new double[nocca*nvirt]);
  fill(lia.get(), lia.get()+nocca*nvirt, 0.0);
  unique_ptr<double[]> lif(new double[nocc*max(static_cast<unsigned long int>(ncore),1lu)]);
  const unique_ptr<double[]> lip = gip->form_2index(geom_->df(), 2.0);
  {
    dgemm_("N", "N", nocc, nvirt, nbasis, 1.0, lip.get(), nocc, vcoeff, nbasis, 0.0, lia.get()+ncore, nocca);
    if (ncore)
      dgemm_("N", "N", nocc, ncore, nbasis, 1.0, lip.get(), nocc, ocoeff, nbasis, 0.0, lif.get(), nocc);
  }

  // core-occ density matrix elements
  for (int i = 0; i != ncore; ++i)
    for (int j = ncore; j != nocca; ++j)
      dmp2->element(j,i) = dmp2->element(i,j) = lif[(j-ncore)+nocc*i] / (eig_tm[j]-eig_tm[i]);

  // 2*J_al(d_rs)
  shared_ptr<Matrix> dmp2ao_part(new Matrix(*ref_->coeff() * *dmp2 ^ *ref_->coeff()));
  unique_ptr<double[]> jai(new double[nvirt*nocca]);
  {
    unique_ptr<double[]> jrs = geom_->df()->compute_Jop(dmp2ao_part->data());
    unique_ptr<double[]> jri(new double[nbasis*nocca]);
    dgemm_("N", "N", nbasis, nocca, nbasis, 1.0, jrs.get(), nbasis, ocoeff, nbasis, 0.0, jri.get(), nbasis);
    dgemm_("T", "N", nvirt, nocca, nbasis, 2.0, vcoeff, nbasis, jri.get(), nbasis, 0.0, jai.get(), nvirt);
  }
  // -1*K_al(d_rs)
  unique_ptr<double[]> kia(new double[nvirt*nocca]);
  {
    unique_ptr<double[]> kir = halfjj->compute_Kop_1occ(dmp2ao_part->data());
    dgemm_("N", "N", nocca, nvirt, nbasis, -1.0, kir.get(), nocca, vcoeff, nbasis, 0.0, kia.get(), nocca);
  }

  shared_ptr<Matrix> grad(new Matrix(nbasis, nbasis));
  for (int i = 0; i != nocca; ++i)
    for (int a = 0; a != nvirt; ++a)
      // minus sign is due to the convention in the solvers which solve Ax+B=0..
      grad->element(a+nocca, i) = - (lai[a+nvirt*i] - lia[i+nocca*a] - jai[a+nvirt*i] - kia[i+nocca*a]);

  auto tp4 = chrono::high_resolution_clock::now();
  auto dr3 = chrono::duration_cast<chrono::milliseconds>(tp4-tp3);
  cout << setw(60) << left << "    * Right hand side of CPHF done" << right << setw(10) << setprecision(2) << dr3.count()*0.001 << endl << endl;

  // solving CPHF
  shared_ptr<CPHF> cphf(new CPHF(grad, ref_->eig(), halfjj, ref_));
  shared_ptr<Matrix> dia = cphf->solve();
  *dmp2 += *dia;

  // total density matrix
  shared_ptr<Matrix> dtot(new Matrix(*dmp2));
  for (int i = 0; i != nocca; ++i) dtot->element(i,i) += 2.0;

  // computes dipole mements
  shared_ptr<Matrix> dtotao(new Matrix(*ref_->coeff() * *dtot ^ *ref_->coeff()));
  Dipole dipole(geom_, dtotao);
  dipole.compute();

  ////////////////////////////////////////////////////////////////////////////

  auto tp5 = chrono::high_resolution_clock::now();
  auto dr4 = chrono::duration_cast<chrono::milliseconds>(tp5-tp4);
  cout << setw(60) << left << "    * CPHF solved" << right << setw(10) << setprecision(2) << dr4.count()*0.001 << endl;

  // one electron matrices
  shared_ptr<Matrix> dmp2ao(new Matrix(*ref_->coeff() * *dmp2 ^ *ref_->coeff()));
  shared_ptr<Matrix> d0ao(new Matrix(*dtotao - *dmp2ao));
  shared_ptr<Matrix> dbarao(new Matrix(*dtotao - *d0ao*0.5));

  // size of naux
  unique_ptr<double[]> cd0 = geom_->df()->compute_cd(d0ao->data());
  unique_ptr<double[]> cdbar = geom_->df()->compute_cd(dbarao->data());


  // three-index derivatives (seperable part)...
  vector<const double*> cd; cd.push_back(cd0.get());      cd.push_back(cdbar.get());
  vector<const double*> dd; dd.push_back(dbarao->data()); dd.push_back(d0ao->data());
  shared_ptr<DF_AO> sep3(new DF_AO(nbasis, nbasis, naux, cd, dd));

  shared_ptr<DF_Half> sepd = halfjj->apply_density(dbarao->data());
  {
    shared_ptr<DF_AO> sep32 = sepd->back_transform(ocoeff);
    sep3->daxpy(-2.0, sep32);
  }
  /// mp2 two body part ----------------
  {
    shared_ptr<DF_AO> sep32 = gip->back_transform(coeff);
    sep3->daxpy(4.0, sep32);
  }

  // two-index derivatives (seperable part)..
  shared_ptr<Matrix> sep2(new Matrix(naux, naux));
  dger_(naux, naux, 2.0, cd0.get(), 1, cdbar.get(), 1, sep2->data(), naux);
  *sep2 -= *halfjj->form_aux_2index(sepd) * 2.0;
  *sep2 += *gia->form_aux_2index_apply_J(full) * 4.0;


  // energy weighted density
  shared_ptr<Matrix> wd(new Matrix(nbasis, nbasis));
  for (int i = 0; i != nocca; ++i)
    for (int j = 0; j != nocca; ++j)
      wd->element(j,i) += dtot->element(j,i) * eig_tm[j];
  for (int i = 0; i != nvirt; ++i)
    for (int j = 0; j != nvirt; ++j)
      wd->element(j+nocca,i+nocca) += dmp2->element(j+nocca,i+nocca) * eig_tm[j+nocca];
  for (int i = 0; i != nocca; ++i)
    for (int j = 0; j != nvirt; ++j)
      wd->element(j+nocca,i) += 2.0 * dmp2->element(j+nocca,i) * eig_tm[i];
  // Liq + Laq
  dgemm_("N", "N", nocc, nocca, nbasis, 1.0, lip.get(), nocc, ocoeff, nbasis, 1.0, wd->data()+ncore, nbasis);
  dgemm_("N", "N", nvirt, nocca, nbasis, 2.0, laq.get(), nvirt, ocoeff, nbasis, 1.0, wd->data()+nocca, nbasis);
  dgemm_("N", "N", nvirt, nvirt, nbasis, 1.0, laq.get(), nvirt, vcoeff, nbasis, 1.0, wd->data()+nocca+nocca*nbasis, nbasis);

  unique_ptr<double[]> jrs = geom_->df()->compute_Jop(dmp2ao->data());
  unique_ptr<double[]> jri(new double[nbasis*nocca]);
  dgemm_("N", "N", nbasis, nocca, nbasis, 1.0, jrs.get(), nbasis, ocoeff, nbasis, 0.0, jri.get(), nbasis);
  dgemm_("T", "N", nocca, nocca, nbasis, 2.0, ocoeff, nbasis, jri.get(), nbasis, 1.0, wd->data(), nbasis);
  unique_ptr<double[]> kir = halfjj->compute_Kop_1occ(dmp2ao->data());
  dgemm_("N", "N", nocca, nocca, nbasis, -1.0, kir.get(), nocca, ocoeff, nbasis, 1.0, wd->data(), nbasis);

  wd->symmetrize();
  shared_ptr<Matrix> wdao(new Matrix(*ref_->coeff() * *wd ^ *ref_->coeff()));

  auto tp6 = chrono::high_resolution_clock::now();
  auto dr5 = chrono::duration_cast<chrono::milliseconds>(tp6-tp5);
  cout << setw(60) << left << "    * Density matrices computed" << right << setw(10) << setprecision(2) << dr5.count()*0.001 << endl;

  // gradient evaluation
  shared_ptr<GradFile> gradf = contract_gradient(dtotao, wdao, sep3, sep2);

  auto tp7 = chrono::high_resolution_clock::now();
  auto dr6 = chrono::duration_cast<chrono::milliseconds>(tp7-tp6);
  cout << setw(60) << left << "    * Gradient integrals contracted " << setprecision(2) << right << setw(10) << dr6.count()*0.001 << endl << endl;

  // set proper energy_
  energy_ = ref_->energy() + ecorr;

  return gradf;

}


