//
// BAGEL - Parallel electron correlation program.
// Filename: mp2.cc
// Copyright (C) 2012 Toru Shiozaki
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


// implements the MP2-F12 theory

#include <stddef.h>
#include <src/mp2/mp2.h>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/scf/scf.h>
#include <src/smith/prim_op.h>
#include <src/mp2/f12int4.h>
#include <src/util/taskqueue.h>
#include <src/parallel/resources.h>

using namespace std;
using namespace bagel;

MP2::MP2(const shared_ptr<const PTree> input, const shared_ptr<const Geometry> g, const shared_ptr<const Reference> ref) : Method(input, g, ref) {

  scf_ = make_shared<SCF>(input, g, ref);
  scf_->compute();
  ref_ = scf_->conv_to_ref();

  cout << endl << "  === DF-MP2 calculation ===" << endl << endl;

  // checks for frozen core
  const bool frozen = idata_->get<bool>("frozen", false);
  ncore_ = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (ncore_) cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;

  if (geom_->df() == nullptr) throw logic_error("MP2 is only implemented with DF");

  // if three is a aux_basis keyword, we use that basis
  abasis_ = to_lower(idata_->get<string>("aux_basis", ""));

}


void MP2::compute() {
  const size_t nbasis = ref_->coeff()->mdim();
  const size_t nocc = ref_->nocc() - ncore_;
  if (nocc < 1) throw runtime_error("no correlated electrons");
  const size_t nvirt = nbasis - nocc - ncore_;
  if (nvirt < 1) throw runtime_error("no virtuals orbitals");


  shared_ptr<const Matrix> ocoeff = ref_->coeff()->slice(ncore_, ncore_+nocc);
  shared_ptr<const Matrix> vcoeff = ref_->coeff()->slice(ncore_+nocc, ncore_+nocc+nvirt);

  Timer timer;

  // first compute half transformed integrals
  shared_ptr<DFHalfDist> half;
  if (abasis_.empty()) {
    half = geom_->df()->compute_half_transform(ocoeff);
  } else {
    auto info = make_shared<PTree>(); info->put("df_basis", abasis_);
    auto cgeom = make_shared<Geometry>(*geom_, info, false);
    half = cgeom->df()->compute_half_transform(ocoeff);
  }
  // second transform for virtual index
  // this is now (naux, nocc, nvirt)
  shared_ptr<DFFullDist> full = half->compute_second_transform(vcoeff)->apply_J();

  cout << "    * 3-index integral transformation done" << endl;

  // assemble
  vector<double> eig(ref_->eig().begin()+ncore_, ref_->eig().end());
  auto buf = make_shared<Matrix>(nocc*nvirt, nocc); // it is implicitly assumed that o^2v can be kept in core in each node
  energy_ = 0.0;
  for (size_t i = 0; i != nvirt; ++i) {
    shared_ptr<const Matrix> data = full->form_4index_1fixed(full, 1.0, i);
    *buf = *data;
    // using SMITH's symmetrizer (src/smith/prim_op.h)
    SMITH::sort_indices<2,1,0,2,1,-1,1>(data->data(), buf->data(), nocc, nvirt, nocc);
    double* tdata = buf->data();
    for (size_t j = 0; j != nocc; ++j)
      for (size_t k = 0; k != nvirt; ++k)
        for (size_t l = 0; l != nocc; ++l, ++tdata)
          *tdata /= -eig[i+nocc]+eig[j]-eig[k+nocc]+eig[l];
    energy_ += data->dot_product(buf);
  }

  cout << "    * assembly done" << endl << endl;
  cout << "      MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << energy_ << setw(10) << setprecision(2) << timer.tick() << endl << endl;

  energy_ += ref_->energy();
  cout << "      MP2 total energy:       " << fixed << setw(15) << setprecision(10) << energy_ << endl << endl;

  // check if F12 is requested.
  const bool do_f12 = idata_->get<bool>("f12", false);
  if (do_f12) {
#ifdef HAVE_LIBSLATER
    const double gamma = idata_->get<double>("gamma", 1.5);
    cout << "    * F12 calculation requested with gamma = " << setprecision(2) << gamma << endl;
#if 0
    auto f12int = make_shared<F12Int>(idata_, geom_, ref_, gamma, ncore_);
#else
    auto f12ref = make_shared<F12Ref>(geom_, ref_, ncore_, gamma);
    f12ref->compute();
#endif
#else
  throw runtime_error("Slater-quadrature library not linked");
#endif
  }
}


