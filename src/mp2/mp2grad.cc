//
// Newint - Parallel electron correlation program.
// Filename: mp2grad.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/mp2/mp2grad.h>
#include <src/grad/cphf.h>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/smith/prim_op.h>
#include <src/prop/dipole.h>
#include <src/grad/gradeval_base.h>
#include <src/util/linear.h>

using namespace std;

MP2Grad::MP2Grad(const multimap<string, string> input, const shared_ptr<Geometry> g, shared_ptr<Reference> r) : MP2(input, g, r) {

}

void MP2Grad::compute() {
  size_t time = ::clock();

  // since this is only for closed shell
  const size_t naux = geom_->naux();
  const size_t nocc = ref_->nocc() - ncore_;
  const size_t nocca = ref_->nocc();
  if (nocc < 1) throw runtime_error("no correlated electrons"); 
  const size_t nvirt = geom_->nbasis() - nocca;
  if (nvirt < 1) throw runtime_error("no virtuals orbitals"); 
  assert(geom_->nbasis() == ref_->coeff()->mdim());

  const size_t nbasis = geom_->nbasis();

  const double* const coeff = ref_->coeff()->data() + ncore_*nbasis;
  const double* const ocoeff = ref_->coeff()->data();
  const double* const vcoeff = coeff + nocc*nbasis;

  // first compute half transformed integrals
  shared_ptr<DF_Half> half = geom_->df()->compute_half_transform(coeff, nocc);  
  // TODO this is a waste...
  shared_ptr<DF_Half> halfjj = geom_->df()->compute_half_transform(ocoeff, nocca)->apply_JJ();
  // second transform for virtual index
  // this is now (naux, nocc, nvirt)
  shared_ptr<DF_Full> full = half->compute_second_transform(vcoeff, nvirt)->apply_J();
  shared_ptr<DF_Full> bv = full->apply_J();
  shared_ptr<DF_Full> gia = bv->clone();

  double elapsed = (::clock()-time)/static_cast<double>(CLOCKS_PER_SEC); 
  cout << setw(60) << left << "    * 3-index integral transformation done" << right << setw(10) << setprecision(2) << elapsed << endl << endl;
  time = ::clock();

  // assemble
  unique_ptr<double[]> buf(new double[nocc*nvirt*nocc]); // it is implicitly assumed that o^2v can be kept in core in each node
  unique_ptr<double[]> buf2(new double[nocc*nvirt*nocc]);
  vector<double> eig_tm = ref_->eig();
  vector<double> eig(eig_tm.begin()+ncore_, eig_tm.end());

  shared_ptr<Matrix1e> dmp2(new Matrix1e(geom_));
  double* optr = dmp2->element_ptr(ncore_, ncore_);
  double* vptr = dmp2->element_ptr(nocca, nocca);

  double sum = 0.0;
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
    sum += ddot_(nocc*nvirt*nocc, data, 1, buf, 1);

    // form Gia : TODO distribute
    // Gia(D|ic) = BV(D|ja) G_c(ja|i)
    dgemm_("N", "N", naux, nocc, nvirt*nocc, 1.0, bv->data(), naux, buf.get(), nocc*nvirt, 0.0, gia->data()+i*nocc*naux, naux); 

    // G(ja|ic) -> G_c(a,ij) 
    SMITH::sort_indices<1,2,0,0,1,1,1>(buf, data, nocc, nvirt, nocc);
    // T(jb|ic) -> T_c(b,ij) 
    SMITH::sort_indices<1,2,0,0,1,1,1>(buf2, buf, nocc, nvirt, nocc);
    // D_ab = G(ja|ic) T(jb|ic)
    dgemm_("N", "T", nvirt, nvirt, nocc*nocc, 2.0, buf.get(), nvirt, data.get(), nvirt, 1.0, vptr, nbasis); 
    // D_ij = - G(ja|kc) T(ia|kc)
    dgemm_("T", "N", nocc, nocc, nvirt*nocc, -2.0, buf.get(), nvirt*nocc, data.get(), nvirt*nocc, 1.0, optr, nbasis); 
  }

  elapsed = (::clock()-time)/static_cast<double>(CLOCKS_PER_SEC); 
  cout << left << setw(60) << "    * assembly (+ unrelaxed density matrices) done" << right << setw(10) << setprecision(2) << elapsed << endl << endl;
  cout << "      MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << sum << endl << endl;
  time = ::clock();

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
  unique_ptr<double[]> lif(new double[nocc*max(ncore_,1)]);
  const unique_ptr<double[]> lip = gip->form_2index(geom_->df(), 2.0);
  {
    dgemm_("N", "N", nocc, nvirt, nbasis, 1.0, lip.get(), nocc, vcoeff, nbasis, 0.0, lia.get()+ncore_, nocca); 
    if (ncore_)
      dgemm_("N", "N", nocc, ncore_, nbasis, 1.0, lip.get(), nocc, ocoeff, nbasis, 0.0, lif.get(), nocc); 
  }

  // core-occ density matrix elements
  for (int i = 0; i != ncore_; ++i)
    for (int j = ncore_; j != nocca; ++j)
      dmp2->element(j,i) = dmp2->element(i,j) = lif[(j-ncore_)+nocc*i] / (eig_tm[j]-eig_tm[i]);

  // 2*J_al(d_rs)
  shared_ptr<Matrix1e> dmp2ao_part(new Matrix1e(*ref_->coeff() * *dmp2 ^ *ref_->coeff()));
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

  shared_ptr<Matrix1e> grad(new Matrix1e(geom_));
  for (int i = 0; i != nocca; ++i)
    for (int a = 0; a != nvirt; ++a)
      // minus sign is due to the convention in the solvers which solve Ax+B=0..
      grad->element(a+nocca, i) = - (lai[a+nvirt*i] - lia[i+nocca*a] - jai[a+nvirt*i] - kia[i+nocca*a]);

  elapsed = (::clock()-time)/static_cast<double>(CLOCKS_PER_SEC); 
  cout << setw(60) << left << "    * Right hand side of CPHF done" << right << setw(10) << setprecision(2) << elapsed << endl << endl;
  time = ::clock();

  // solving CPHF
  shared_ptr<CPHF> cphf(new CPHF(grad, ref_->eig(), halfjj, ref_));
  shared_ptr<Matrix1e> dia = cphf->solve();
  *dmp2 += *dia;

  // total density matrix
  shared_ptr<Matrix1e> dtot(new Matrix1e(*dmp2));
  for (int i = 0; i != nocca; ++i) dtot->element(i,i) += 2.0;

  // computes dipole mements
  shared_ptr<Matrix1e> dtotao(new Matrix1e(*ref_->coeff() * *dtot ^ *ref_->coeff()));
  Dipole dipole(geom_, dtotao);
  dipole.compute();

// dipole is correct
////////////////////////////////////////////////////////////////////////////

  elapsed = (::clock()-time)/static_cast<double>(CLOCKS_PER_SEC); 
  cout << endl;
  cout << setw(60) << left << "    * CPHF solved" << right << setw(10) << setprecision(2) << elapsed << endl << endl;
  time = ::clock();

  // one electron matrices
  shared_ptr<Matrix1e> dmp2ao(new Matrix1e(*ref_->coeff() * *dmp2 ^ *ref_->coeff()));
  shared_ptr<Matrix1e> d0ao(new Matrix1e(*dtotao - *dmp2ao));
  shared_ptr<Matrix1e> dbarao(new Matrix1e(*dtotao - *d0ao*0.5));

  // size of naux
  unique_ptr<double[]> cd0 = geom_->df()->compute_cd(d0ao->data()); 
  unique_ptr<double[]> cdbar = geom_->df()->compute_cd(dbarao->data()); 


  // three-index derivatives (seperable part)...
  vector<const double*> cd; cd.push_back(cd0.get());      cd.push_back(cdbar.get());
  vector<const double*> dd; dd.push_back(dbarao->data()); dd.push_back(d0ao->data()); 
  shared_ptr<DF_AO> sep3(new DF_AO(nbasis, nbasis, naux, cd, dd));

  // TODO perhaps we could merge back transformation...
  shared_ptr<DF_Half> sepd = halfjj->apply_density(dbarao->data());
  {
    shared_ptr<DF_AO> sep32 = sepd->back_transform(ref_->coeff()->data());
    sep3->daxpy(-2.0, sep32);
  }
  {
    shared_ptr<DF_AO> sep32 = gip->back_transform(ref_->coeff()->data());
    sep3->daxpy(1.0, sep32);
  }
  // symmetrize...
  for (int i = 0; i != nbasis; ++i) {
    for (int j = i+1; j != nbasis; ++j) {
      for (int k = 0; k != naux; ++k) {
        const double tmp = 0.5*(*sep3->ptr(k,j,i) + *sep3->ptr(k,i,j));
        *sep3->ptr(k,j,i) = *sep3->ptr(k,i,j) = tmp;
      }
    }
  }

  // two-index derivatives (seperable part)..
  unique_ptr<double[]> sep2(new double[naux*naux]);
  fill(sep2.get(), sep2.get()+naux*naux, 0.0);
  dger_(naux, naux, 2.0, cd0, 1, cdbar, 1, sep2, naux); 
  {
    unique_ptr<double[]> sep22 = halfjj->form_aux_2index(sepd);
    daxpy_(naux*naux, -2.0, sep22, 1, sep2, 1);
  }
  {
    unique_ptr<double[]> sep22 = gia->form_aux_2index(full);
    dgemm_("N", "N", naux, naux, naux, -2.0, sep22.get(), naux, geom_->df()->data_2index(), naux, 1.0, sep2.get(), naux);
  }

  // symmetrize..
  for (int i = 0; i != naux; ++i)
    for (int j = i+1; j != naux; ++j)
      sep2[j+i*naux] = sep2[i+j*naux] = 0.5*(sep2[j+i*naux] + sep2[i+j*naux]); 


  // energy weighted density
  shared_ptr<Matrix1e> wd(new Matrix1e(geom_));
  for (int i = 0; i != nocc; ++i)
    for (int j = 0; j != nocc; ++j)
      wd->element(j,i) += 0.5 * dtot->element(j,i) * (eig[j] + eig[i]); 
  for (int i = 0; i != nvirt; ++i)
    for (int j = 0; j != nvirt; ++j)
      wd->element(j+nocc,i+nocc) += 0.5 * dtot->element(j+nocc,i+nocc) * (eig[j+nocc] + eig[i+nocc]); 
  for (int i = 0; i != nocc; ++i)
    for (int j = 0; j != nvirt; ++j)
      wd->element(j+nocc,i) += 2.0*dtot->element(j+nocc,i) * eig[i];
  // Liq + Laq
  dgemm_("N", "N", nocc, nocc, nbasis, 1.0, lip.get(), nocc, ocoeff, nbasis, 1.0, wd->data(), nbasis); 
  dgemm_("N", "N", nvirt, nocc, nbasis, 2.0, laq.get(), nvirt, ocoeff, nbasis, 1.0, wd->data()+nocc, nbasis); 
  dgemm_("N", "N", nvirt, nvirt, nbasis, 1.0, laq.get(), nvirt, vcoeff, nbasis, 1.0, wd->data()+nocc+nocc*nbasis, nbasis); 

  unique_ptr<double[]> jrs = geom_->df()->compute_Jop(dmp2ao->data());
  unique_ptr<double[]> jri(new double[nbasis*nocc]);
  dgemm_("N", "N", nbasis, nocc, nbasis, 1.0, jrs.get(), nbasis, coeff, nbasis, 0.0, jri.get(), nbasis);
  dgemm_("T", "N", nocc, nocc, nbasis, 2.0, coeff, nbasis, jri.get(), nbasis, 1.0, wd->data(), nbasis); 
  unique_ptr<double[]> kir = halfjj->compute_Kop_1occ(dmp2ao->data());
  dgemm_("N", "N", nocc, nocc, nbasis, -1.0, kir.get(), nocc, coeff, nbasis, 1.0, wd->data(), nbasis); 

  wd->symmetrize();
  wd->print();
  shared_ptr<Matrix1e> wdao(new Matrix1e(*ref_->coeff() * *wd ^ *ref_->coeff()));

  // gradient evaluation
  GradEval_base eval(geom_);
  eval.contract_gradient(dtotao, wdao, sep3, sep2)->print(); 

}


