//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mp2grad.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/pt2/mp2/mp2grad.h>
#include <src/pt2/mp2/mp2cache.h>
#include <src/pt2/mp2/mp2accum.h>
#include <src/grad/cphf.h>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/util/prim_op.h>
#include <src/prop/multipole.h>
#include <src/grad/gradeval.h>
#include <src/prop/moprint.h>

using namespace std;
using namespace bagel;
using namespace btas;

MP2Grad::MP2Grad(shared_ptr<const PTree> input, shared_ptr<const Geometry> g, shared_ptr<const Reference> ref) : MP2(input, g, ref) {

}


// does nothing.
void MP2Grad::compute() { }


template<>
vector<double> GradEval<MP2Grad>::energyvec() const {
  vector<double> out;
  out.push_back(energy());
  return out;
}


template<>
shared_ptr<GradFile> GradEval<MP2Grad>::compute(const string jobtitle, shared_ptr<const GradInfo> gradinfo) {
  Timer time;

  const size_t ncore = task_->ncore();
  const size_t nmobasis = ref_->coeff()->mdim();

  // since this is only for closed shell
  const size_t nocc = ref_->nocc() - ncore;
  const size_t nocca = ref_->nocc();
  if (nocc < 1) throw runtime_error("no correlated electrons");
  const size_t nvirt = nmobasis - nocca;
  if (nvirt < 1) throw runtime_error("no virtuals orbitals");

  const MatView ccmat = ref_->coeff()->slice(0, ncore);
  const MatView ocmat = ref_->coeff()->slice(0, nocca);
  const MatView vcmat = ref_->coeff()->slice(nocca, nmobasis);
  const MatView acmat = ref_->coeff()->slice(ncore, nocca);

  // first compute half transformed integrals
  shared_ptr<DFHalfDist> half;
  shared_ptr<const Geometry> cgeom;
  if (task_->abasis().empty()) {
    cgeom = geom_;
    half = geom_->df()->compute_half_transform(acmat);
  } else {
    auto info = make_shared<PTree>(); info->put("df_basis", task_->abasis());
    cgeom = make_shared<Geometry>(*geom_, info, false);
    half = cgeom->df()->compute_half_transform(acmat);
  }
  shared_ptr<const DFHalfDist> halfjj = geom_->df()->compute_half_transform(ocmat)->apply_JJ();
  // second transform for virtual index
  // this is now (naux, nocc, nvirt)
  shared_ptr<const DFFullDist> full = half->compute_second_transform(vcmat)->apply_J();
  shared_ptr<const DFFullDist> bv = full->apply_J();
  shared_ptr<DFFullDist> gia = bv->clone();

  time.tick_print("3-index integral transform");

  // assemble
  const VectorB eig_tm = ref_->eig();
  const double* eig = eig_tm.data() + ncore;

  auto dmp2 = make_shared<Matrix>(nmobasis, nmobasis);
  double ecorr = 0.0;
  {
    // second transform for virtual index and rearrange data
    auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size(), full->nocc2());
    auto fullt = make_shared<DFDistT>(full->swap(), dist);

    MP2Cache cache(fullt->naux(), nocc, nvirt, fullt);

    size_t memory_size = half->block(0)->size() * 2;
    mpi__->broadcast(&memory_size, 1, 0);
    const int nloop = cache.nloop();
    const int ncache = min(memory_size/(nvirt*nvirt), size_t(20));
    cout << "    * ncache = " << ncache << endl;
    for (int n = 0; n != min(ncache, nloop); ++n)
      cache.block(n, -1);

    for (int n = 0; n != nloop; ++n) {
      // take care of data. The communication should be hidden
      if (n+ncache < nloop)
        cache.block(n+ncache, n-1);

      const int i = get<0>(cache.task(n));
      const int j = get<1>(cache.task(n));
      if (i < 0 || j < 0) continue;
      cache.data_wait(n);

      shared_ptr<const Matrix> iblock = cache(i);
      shared_ptr<const Matrix> jblock = cache(j);
      const Matrix mat(*iblock % *jblock); // V
      Matrix mat2 = mat; // 2T-T^t
      if (i != j) {
        mat2 *= 2.0;
        mat2 -= *mat.transpose();
      }
      Matrix mat3 = mat; // T

      for (int a = 0; a != nvirt; ++a) {
        for (int b = 0; b != nvirt; ++b) {
          const double denom = -eig[a+nocc]+eig[i]-eig[b+nocc]+eig[j];
          mat2(b,a) /= denom;
          mat3(b,a) /= denom;
        }
      }

      const double fac = i == j ? 1.0 : 2.0;
      ecorr += mat.dot_product(mat2) * fac;
      dmp2->add_block(2.0, nocca, nocca, nvirt, nvirt, mat2 % mat3);
      if (i != j)
        dmp2->add_block(2.0, nocca, nocca, nvirt, nvirt, mat2 ^ mat3);
    }
    cache.wait();
    // allreduce energy contributions
    mpi__->allreduce(&ecorr, 1);
  }

  time.tick_print("First pass based on occupied orbitals");

  // do the same with virtual-virtual index distribution to obtain D_ij
  // TODO is there any better way??
  {
    // second transform for virtual index and rearrange data
    auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size(), full->nocc1());
    auto fullt = make_shared<DFDistT>(full, dist);
    MP2Cache cache(fullt->naux(), nvirt, nocc, fullt);

    auto giat = fullt->clone();
    MP2Accum accum(fullt->naux(), nvirt, nocc, giat, cache.tasks());

    size_t memory_size = half->block(0)->size() * 2;
    mpi__->broadcast(&memory_size, 1, 0);
    const int nloop = cache.nloop();
    const int ncache = min(memory_size/(nocc*nocc), size_t(30));
    cout << "    * ncache = " << ncache << endl;
    for (int n = 0; n != min(ncache, nloop); ++n)
      cache.block(n, -1);

    for (int n = 0; n != nloop; ++n) {
      // take care of data. The communication should be hidden
      if (n+ncache < nloop)
        cache.block(n+ncache, n-1);

      const int i = get<0>(cache.task(n));
      const int j = get<1>(cache.task(n));
      if (i >= 0 && j >= 0) {
        cache.data_wait(n);

        shared_ptr<const Matrix> iblock = cache(i);
        shared_ptr<const Matrix> jblock = cache(j);
        Matrix mat2 = *iblock % *jblock; // 2T-T^t
        Matrix mat3 = mat2; // T
        if (i != j) {
          mat2 *= 2.0;
          mat2 -= *mat3.transpose();
        }

        for (int a = 0; a != nocc; ++a) {
          for (int b = 0; b != nocc; ++b) {
            const double denom = -eig[a]+eig[i+nocc]-eig[b]+eig[j+nocc];
            mat2(b,a) /= denom;
            mat3(b,a) /= denom;
          }
        }

        accum.accumulate<0>(n, make_shared<Matrix>(*jblock ^ mat2));
        accum.accumulate<1>(n, make_shared<Matrix>(*iblock * mat2));

        dmp2->add_block(-2.0, ncore, ncore, nocc, nocc, mat2 % mat3);
        if (i != j)
          dmp2->add_block(-2.0, ncore, ncore, nocc, nocc, mat2 ^ mat3);
      } else {
        // for receiving data
        accum.accumulate<0>(n, nullptr);
        accum.accumulate<1>(n, nullptr);
      }
    }
    accum.wait();
    cache.wait();
    giat = giat->apply_J();
    giat->get_paralleldf(gia);
  }

  dmp2->allreduce();
  gia->scale(-1.0);

  time.tick_print("Second pass based on virtual orbitals");
  cout << endl;
  cout << "      MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << ecorr << endl << endl;

  // L''aq = 2 Gia(D|ia) (D|iq)
  shared_ptr<const Matrix> laq = gia->form_2index(half, 2.0);
  const Matrix lai = *laq * ocmat;

  // Gip = Gia(D|ia) C+_ap
  shared_ptr<DFHalfDist> gip = gia->back_transform(vcmat);
  // Liq = 2 Gip(D) * (D|pq)
  Matrix lia(nocca, nvirt);
  Matrix lif(nocc, max(1lu,ncore));
  shared_ptr<const Matrix> lip = gip->form_2index(cgeom->df(), 2.0);
  {
    lia.add_block(1.0, ncore, 0, nocc, nvirt, *lip * vcmat);
    if (ncore)
      lif = *lip * ccmat;
  }

  // core-occ density matrix elements
  for (int i = 0; i != ncore; ++i)
    for (int j = ncore; j != nocca; ++j)
      dmp2->element(j,i) = dmp2->element(i,j) = lif(j-ncore, i) / (eig_tm(j)-eig_tm(i));

  // 2*J_al(d_rs)
  auto dmp2ao_part = make_shared<const Matrix>(*ref_->coeff() * *dmp2 ^ *ref_->coeff());
  const Matrix jai = vcmat % *geom_->df()->compute_Jop(dmp2ao_part) * ocmat * 2.0;
  // -1*K_al(d_rs)
  const Matrix kia = *halfjj->compute_Kop_1occ(dmp2ao_part, -1.0) * vcmat;

  auto grad = make_shared<Matrix>(nmobasis, nmobasis);
  for (int i = 0; i != nocca; ++i)
    for (int a = 0; a != nvirt; ++a)
      // minus sign is due to the convention in the solvers which solve Ax+B=0..
      grad->element(a+nocca, i) = - (lai(a,i) - lia(i,a) - jai(a,i) - kia(i,a));

  {
    auto d_unrelaxed = make_shared<Matrix>(*dmp2);
    for (int i = 0; i != nocca; ++i) d_unrelaxed->element(i,i) += 2.0;
    auto dao_unrelaxed = make_shared<Matrix>(*ref_->coeff() * *d_unrelaxed ^ *ref_->coeff());
    Dipole dipole(geom_, dao_unrelaxed, "MP2 unrelaxed");
    dipole.compute();
  }
  time.tick_print("Right hand side of CPHF");
  cout << endl;

  // solving CPHF (or Z-vector equation)
  auto cphf = make_shared<CPHF>(grad, ref_->eig(), halfjj, ref_);
  shared_ptr<Matrix> dia = cphf->solve(task_->scf()->thresh_scf(), gradinfo->maxziter());
  *dmp2 += *dia;

  // total density matrix
  auto dtot = make_shared<Matrix>(*dmp2);
  for (int i = 0; i != nocca; ++i) dtot->element(i,i) += 2.0;

  // computes dipole mements
  auto dtotao = make_shared<Matrix>(*ref_->coeff() * *dtot ^ *ref_->coeff());
  {
    Dipole dipole(geom_, dtotao, "MP2 relaxed");
    dipole_ = dipole.compute();
  }

  // print relaxed density if requested
  if (gradinfo->density_print()) {
    auto density_print = make_shared<MOPrint>(gradinfo->moprint_info(), geom_, ref_, /*is_density=*/true, make_shared<const ZMatrix>(*dtotao, 1.0));
    density_print->compute();
  }

  ////////////////////////////////////////////////////////////////////////////

  time.tick_print("CPHF solved");

  // one electron matrices
  auto dmp2ao = make_shared<Matrix>(*ref_->coeff() * *dmp2 ^ *ref_->coeff());
  auto d0ao = make_shared<Matrix>(*dtotao - *dmp2ao);
  auto dbarao = make_shared<Matrix>(*dtotao - *d0ao*0.5);

  // size of naux
  shared_ptr<const VectorB> cd0 = geom_->df()->compute_cd(d0ao);
  shared_ptr<const VectorB> cdbar = geom_->df()->compute_cd(dbarao);


  // three-index derivatives (seperable part)...
  vector<shared_ptr<const VectorB>> cd {cd0, cdbar};
  vector<shared_ptr<const Matrix>> dd {dbarao, d0ao};

  shared_ptr<DFHalfDist> sepd = halfjj->apply_density(dbarao);
  shared_ptr<DFDist> sep3 = sepd->back_transform(ocmat);
  sep3->scale(-2.0);
  sep3->add_direct_product(cd, dd, 1.0);
  /// mp2 two body part ----------------
  shared_ptr<DFDist> sep32;
  if (geom_ == cgeom) {
    sep3->ax_plus_y(4.0, gip->back_transform(acmat));
  } else {
    sep32 = gip->back_transform(acmat);
    sep32->scale(4.0);
  }

  // two-index derivatives (seperable part)..
  auto sep2 = make_shared<Matrix>((*cd0 ^ *cdbar) * 2.0);
  *sep2 -= *halfjj->form_aux_2index(sepd, 2.0);
  shared_ptr<const Matrix> sep22;
  if (geom_ == cgeom) {
    *sep2 += *gia->form_aux_2index_apply_J(full, 4.0);
  } else {
    sep22 = gia->form_aux_2index_apply_J(full, 4.0);
  }

  // energy weighted density
  auto wd = make_shared<Matrix>(nmobasis, nmobasis);
  for (int i = 0; i != nocca; ++i)
    for (int j = 0; j != nocca; ++j)
      wd->element(j,i) += dtot->element(j,i) * eig_tm(j);
  for (int i = 0; i != nvirt; ++i)
    for (int j = 0; j != nvirt; ++j)
      wd->element(j+nocca,i+nocca) += dmp2->element(j+nocca,i+nocca) * eig_tm(j+nocca);
  for (int i = 0; i != nocca; ++i)
    for (int j = 0; j != nvirt; ++j)
      wd->element(j+nocca,i) += 2.0 * dmp2->element(j+nocca,i) * eig_tm(i);
  // Liq + Laq
  wd->add_block(1.0, ncore, 0, nocc, nocca, *lip * ocmat);
  wd->add_block(1.0, nocca, 0, nvirt, nocca, *laq * ocmat * 2.0);
  wd->add_block(1.0, nocca, nocca, nvirt, nvirt, *laq * vcmat);
  wd->add_block(1.0, 0, 0, nocca, nocca, (ocmat % *geom_->df()->compute_Jop(dmp2ao) * ocmat * 2.0));
  wd->add_block(1.0, 0, 0, nocca, nocca, (*halfjj->compute_Kop_1occ(dmp2ao, -1.0) * ocmat));

  wd->symmetrize();
  auto wdao = make_shared<Matrix>(*ref_->coeff() * *wd ^ *ref_->coeff());

  time.tick_print("Density matrices computed");
  cout << endl;

  // gradient evaluation
  shared_ptr<GradFile> gradf = contract_gradient(dtotao, wdao, sep3, sep2, /*nacme = */nullptr, false, cgeom, sep32, sep22);

  time.tick_print("Gradient integrals contracted");
  cout << endl;

  // set proper energy_
  energy_ = ref_->energy(0) + ecorr;

  gradf->print();
  return gradf;

}


