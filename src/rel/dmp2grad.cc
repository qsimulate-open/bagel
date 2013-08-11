//
// BAGEL - Parallel electron correlation program.
// Filename: dmp2grad.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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


// implements the 4-component MP2 theory

#include <stddef.h>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/rel/dmp2.h>
#include <src/rel/dfock.h>
#include <src/rel/reldffull.h>
#include <src/rel/relreference.h>
#include <src/rel/reldipole.h>
#include <src/grad/gradeval.h>
#include <src/smith/prim_op.h>
#include <src/parallel/resources.h>

using namespace std;
using namespace bagel;

DMP2Grad::DMP2Grad(const shared_ptr<const PTree> input, const shared_ptr<const Geometry> g, const shared_ptr<const Reference> ref) : DMP2(input, g, ref) {

}


// does nothing.
void DMP2Grad::compute() { }


template<>
shared_ptr<GradFile> GradEval<DMP2Grad>::compute() {
  Timer timer;

  const size_t nbasis = geom_->nbasis();
  shared_ptr<const RelReference> ref = dynamic_pointer_cast<const RelReference>(ref_);

  const size_t ncore = task_->ncore();
  const size_t nocca = ref_->nocc();
  const size_t nocc = nocca - ncore;
  if (nocc < 1) throw runtime_error("no correlated electrons");
  const size_t nvirt = nbasis*2 - nocc - ncore;
  if (nvirt < 1) throw runtime_error("no virtuals orbitals");

  assert(nbasis*4 == ref->relcoeff()->ndim());
  assert(nbasis*2 == ref->relcoeff()->mdim());

  // Separate Coefficients into real and imaginary
  // correlated occupied orbitals
  array<shared_ptr<const Matrix>, 4> rocoeff;
  array<shared_ptr<const Matrix>, 4> iocoeff;
  // correlated virtual orbitals
  array<shared_ptr<const Matrix>, 4> rvcoeff;
  array<shared_ptr<const Matrix>, 4> ivcoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> oc = ref->relcoeff()->get_submatrix(i*nbasis, ncore, nbasis, nocc);
    rocoeff[i] = oc->get_real_part();
    iocoeff[i] = oc->get_imag_part();
    shared_ptr<const ZMatrix> vc = ref->relcoeff()->get_submatrix(i*nbasis, nocca, nbasis, nvirt);
    rvcoeff[i] = vc->get_real_part();
    ivcoeff[i] = vc->get_imag_part();
  }

  // (1) make DFDists
  shared_ptr<Geometry> cgeom;
  vector<shared_ptr<const DFDist>> dfs;
  if (task_->abasis().empty()) {
    dfs = geom_->dfs()->split_blocks();
    dfs.push_back(geom_->df());
  } else {
    auto info = make_shared<PTree>(); info->put("df_basis", task_->abasis());
    cgeom = make_shared<Geometry>(*geom_, info, false);
    cgeom->relativistic(false /* do_gaunt */);
    dfs = cgeom->dfs()->split_blocks();
    dfs.push_back(cgeom->df());
  }
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) first-transform
  list<shared_ptr<RelDFHalf>> half_complex = DFock::make_half_complex(dfdists, rocoeff, iocoeff);

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_exch;
  for (auto& i : half_complex) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
  }
  half_complex.clear();
  DFock::factorize(half_complex_exch);

  // (4) compute (gamma|ia)
  list<shared_ptr<RelDFFull>> dffull;
  for (auto& i : half_complex_exch)
    dffull.push_back(make_shared<RelDFFull>(i, rvcoeff, ivcoeff));
  DFock::factorize(dffull);
  dffull.front()->scale(dffull.front()->fac()); // take care of the factor
  assert(dffull.size() == 1);
  shared_ptr<const RelDFFull> full = dffull.front(); 

  cout << "    * 3-index integral transformation done" << endl;

  // assemble
  vector<double> eig(ref_->eig().begin()+ncore, ref_->eig().end());
  auto buf = make_shared<ZMatrix>(nocc*nvirt, nocc); // it is implicitly assumed that o^2v can be kept in core in each node

  // (\check{gamam}|ia)
  shared_ptr<RelDFFull> bv = full->apply_J();
  shared_ptr<RelDFFull> gia = bv->clone();

  auto dmp2 = make_shared<ZMatrix>(nbasis*2, nbasis*2);
  complex<double>* optr = dmp2->element_ptr(ncore, ncore);
  complex<double>* vptr = dmp2->element_ptr(nocca, nocca);

  energy_ = 0.0;
  for (size_t i = 0; i != nvirt; ++i) {
    shared_ptr<ZMatrix> data = full->form_4index_1fixed(full, 1.0, i);
    *buf = *data;
    // using SMITH's symmetrizer (src/smith/prim_op.h)
    SMITH::sort_indices<2,1,0,1,1,-1,1>(data->data(), buf->data(), nocc, nvirt, nocc);
    complex<double>* tdata = buf->data();
    for (size_t j = 0; j != nocc; ++j)
      for (size_t k = 0; k != nvirt; ++k)
        for (size_t l = 0; l != nocc; ++l, ++tdata)
          *tdata /= -eig[i+nocc]+eig[j]-eig[k+nocc]+eig[l]; // assumed that the denominator is positive
    energy_ += data->dot_product(buf).real() * 0.5;

    // form Gia : TODO distribute
    // Gia(D|ic) = BV(D|ja) G_c(ja|i)
    // BV and gia are DFFullDist
#if 0
    const size_t offset = i*nocc;
    gia->add_product(bv, buf, nocc, offset);
#endif

    // T(jb|ic) -> T_c(b,ij)
    SMITH::sort_indices<1,2,0,0,1,1,1>(buf->data(), data->data(), nocc, nvirt, nocc);
    // D_ab = T^*(ja|ic) T(jb|ic)
    zgemm3m_("N", "C", nvirt, nvirt, nocc*nocc, 0.5, data->data(), nvirt, data->data(), nvirt, 1.0, vptr, nbasis*2);
    // D_ij = - T^*(ja|kc) T(ia|kc)
//  std::for_each(data->data(), data->data()+data->size(), [](std::complex<double> a){ a = std::conj(a); });
    zgemm3m_("C", "N", nocc, nocc, nvirt*nocc, -0.5, data->data(), nvirt*nocc, data->data(), nvirt*nocc, 1.0, optr, nbasis*2);
  }

  timer.tick_print("assembly (+ unrelaxed rdm)");
  cout << endl;
  cout << "      DMP2 correlation energy: " << fixed << setw(15) << setprecision(10) << energy_ << setw(10) << setprecision(2) << timer.tick() << endl << endl;

  energy_ += ref_->energy();
  cout << "      DMP2 total energy:       " << fixed << setw(15) << setprecision(10) << energy_ << endl << endl;

  {
    auto d_unrelaxed = make_shared<ZMatrix>(*dmp2);
    for (int i = 0; i != nocca; ++i) d_unrelaxed->element(i,i) += 1.0;
    auto dao_unrelaxed = make_shared<const ZMatrix>(*ref->relcoeff() * *d_unrelaxed ^ *ref->relcoeff());
    RelDipole dipole(geom_, dao_unrelaxed, "MP2 unrelaxed");
    dipole.compute();
  }

  throw logic_error("not yet implemented");
  return shared_ptr<GradFile>();
}


