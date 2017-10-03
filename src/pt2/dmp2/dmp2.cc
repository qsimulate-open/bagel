//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dmp2.cc
// Copyright (C) 2013 Toru Shiozaki
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


// implements the 4-component MP2 theory

#include <stddef.h>
#include <iostream>
#include <iomanip>
#include <src/df/reldffull.h>
#include <src/pt2/dmp2/dmp2.h>
#include <src/scf/dhf/dirac.h>
#include <src/scf/dhf/dfock.h>
#include <src/ci/zfci/reljop.h>
#include <src/wfn/relreference.h>
#include <src/util/prim_op.h>
#include <src/util/f77.h>
#include <src/util/parallel/resources.h>

using namespace std;
using namespace bagel;

DMP2::DMP2(shared_ptr<const PTree> input, shared_ptr<const Geometry> g, shared_ptr<const Reference> ref) : Method(input, g, ref) {

  if (geom_->dfs() || (ref && ref->geom()->dfs())) {
    if (!geom_) geom_ = ref->geom();
  } else {
    auto scf = make_shared<Dirac>(input, g, ref);
    scf->compute();
    ref_ = scf->conv_to_ref();
    geom_ = ref_->geom();
  }

  auto relref = dynamic_pointer_cast<const RelReference>(ref_);
  gaunt_ = relref->gaunt();
  breit_ = relref->breit();

  geom_ = geom_->relativistic(gaunt_);
  assert(geom_->dfs());

  cout << endl << "  === Four-Component DF-MP2 calculation ===" << endl << endl;

  // checks for frozen core
  const bool frozen = idata_->get<bool>("frozen", true);
  ncore_ = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only() : 0));
  if (ncore_) cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;

  // if three is a aux_basis keyword, we use that basis
  abasis_ = to_lower(idata_->get<string>("aux_basis", ""));

}


void DMP2::compute() {
  Timer timer;

  const size_t nbasis = geom_->nbasis();
  shared_ptr<const RelReference> ref = dynamic_pointer_cast<const RelReference>(ref_);

  const size_t nocc = ref_->nocc() - ncore_;
  const size_t nvirt = nbasis*2 - nocc - ncore_;
  if (nocc < 1)  throw runtime_error("no correlated electrons");
  if (nvirt < 1) throw runtime_error("no virtuals orbitals");

  shared_ptr<const ZMatrix> coeff = ref->relcoeff();
  assert(nbasis*4 == coeff->ndim());
  assert(nbasis*2 == coeff->mdim());

  shared_ptr<const Geometry> cgeom = geom_;
  if (!abasis_.empty()) {
    auto info = make_shared<PTree>(); info->put("df_basis", abasis_);
    auto tmp = make_shared<Geometry>(*geom_, info, false);
    tmp->relativistic(gaunt_);
    cgeom = tmp;
  }

  shared_ptr<const ListRelDFFull> fullc, fullg, fullg2;
  {
    list<shared_ptr<RelDFHalf>> half_coulomb;
    tie(half_coulomb, ignore) = RelJop::compute_half(cgeom, coeff->slice_copy(ncore_, ncore_+nocc), false, false);
    fullc = RelJop::compute_full(coeff->slice_copy(ncore_+nocc, ncore_+nocc+nvirt), half_coulomb, true);
  }
  if (gaunt_) {
    list<shared_ptr<RelDFHalf>> half_gaunt, half_gaunt2;
    tie(half_gaunt, half_gaunt2) = RelJop::compute_half(cgeom, coeff->slice_copy(ncore_, ncore_+nocc), gaunt_, breit_);
    fullg = RelJop::compute_full(coeff->slice_copy(ncore_+nocc, ncore_+nocc+nvirt), half_gaunt, true);
    fullg2 = !breit_ ? fullg : RelJop::compute_full(coeff->slice_copy(ncore_+nocc, ncore_+nocc+nvirt), half_gaunt2, false);
  }

  cout << "    * 3-index integral transformation done" << endl;

  // assemble
  vector<double> eig(ref_->eig().begin()+ncore_, ref_->eig().end());
  auto buf = make_shared<ZMatrix>(nocc*nvirt, nocc); // it is implicitly assumed that o^2v can be kept in core in each node

  energy_ = 0.0;
  for (size_t i = 0; i != nvirt; ++i) {
    shared_ptr<ZMatrix> data = fullc->form_4index_1fixed(fullc, 1.0, i);
    if (gaunt_) {
      *data += *fullg->form_4index_1fixed(fullg2, (breit_ ? -0.25 : -1.0), i);
      if (breit_)
        *data += *fullg2->form_4index_1fixed(fullg, (breit_ ? -0.25 : -1.0), i);
    }
    *buf = *data;
    // using a symmetrizer (src/util/prim_op.h)
    sort_indices<2,1,0,1,1,-1,1>(data->data(), buf->data(), nocc, nvirt, nocc);
    complex<double>* tdata = buf->data();
    for (size_t j = 0; j != nocc; ++j)
      for (size_t k = 0; k != nvirt; ++k)
        for (size_t l = 0; l != nocc; ++l, ++tdata)
          *tdata /= -eig[i+nocc]+eig[j]-eig[k+nocc]+eig[l]; // assumed that the denominator is positive
    energy_ += data->dot_product(buf).real() * 0.5;
  }

  cout << "    * assembly done" << endl << endl;
  cout << "      Dirac MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << energy_ << setw(10) << setprecision(2) << timer.tick() << endl << endl;

  energy_ += ref_->energy(0);
  cout << "      Dirac MP2 total energy:       " << fixed << setw(15) << setprecision(10) << energy_ << endl << endl;

}


