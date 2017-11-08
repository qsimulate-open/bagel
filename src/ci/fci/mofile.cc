//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mofile.cc
// Copyright (C) 2011 Toru Shiozaki
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

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <src/ci/fci/mofile.h>
#include <src/scf/hf/fock.h>
#include <src/util/f77.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;


MOFile::MOFile(const shared_ptr<const Reference> ref, const string method)
: geom_(ref->geom()), ref_(ref), coeff_(ref_->coeff()) {

  const bool do_df = geom_->df().get();
  if (!do_df) throw runtime_error("MOFile is implemented only with density fitting");

  hz_ = method == "HZ";
}


MOFile::MOFile(const shared_ptr<const Reference> ref, const shared_ptr<const Matrix> c, const string method)
: hz_(false), geom_(ref->geom()), ref_(ref), coeff_(c) {

  const bool do_df = geom_->df().get();
  if (!do_df) throw runtime_error("MOFile is implemented only with density fitting");

  hz_ = method == "HZ";
}


void MOFile::init(const int nstart, const int nfence, const bool store) {

  // first compute all the AO integrals in core
  nocc_ = nfence - nstart;

  // core energy is set here
  if (nstart != 0) {
    core_fock_ = make_shared<Fock<1>>(geom_, ref_->hcore(), nullptr, coeff_->slice(0,nstart), /*grad*/store, /*rhf*/true);
    shared_ptr<const Matrix> den = coeff_->form_density_rhf(nstart);
    core_energy_ = (*den * (*ref_->hcore()+*core_fock_)).trace() * 0.5;
  } else {
    core_fock_ = ref_->hcore();
    core_energy_ = 0.0;
  }

  // one electron part
  shared_ptr<const Matrix> buf1e = compute_mo1e(nstart, nfence);

  // two electron part.
  // this fills mo2e_1ext_ and returns buf2e which is an ii/ii quantity
  shared_ptr<const Matrix> buf2e = compute_mo2e(nstart, nfence);

  compress_and_set(buf1e, buf2e);
}


void MOFile::compress_and_set(shared_ptr<const Matrix> buf1e, shared_ptr<const Matrix> buf2e) {

  // mo2e is compressed in KH case, not in HZ case
  const int nocc = nocc_;
  sizeij_ = hz_ ? nocc*nocc : nocc*(nocc+1)/2;
  mo2e_ = make_shared<Matrix>(sizeij_, sizeij_, true);

  if (!hz_) {
    int ij = 0;
    for (int i = 0; i != nocc; ++i) {
      for (int j = 0; j <= i; ++j, ++ij) {
        int kl = 0;
        for (int k = 0; k != nocc; ++k)
          for (int l = 0; l <= k; ++l, ++kl)
            mo2e(kl, ij) = buf2e->element(l+k*nocc, j+i*nocc);
      }
    }
  } else {
    // In this case, there is no compression (this is actually necessary)
    // Is currently ordered like (ij|kl), should be ordered like (ik|jl), with the last index moving the fastest
    // Equivalent to             <ik|jl> --> <ij|kl>
    sort_indices<0,2,1,3,0,1,1,1>(buf2e->data(), mo2e_->data(), nocc_, nocc_, nocc_, nocc_);
  }

  // h'kl = hkl - 0.5 sum_j (kj|jl)
  mo1e_ = make_shared<CSymMatrix>(nocc, true);
  int ij = 0;
  for (int i = 0; i != nocc; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      mo1e(ij) = buf1e->element(j,i);
      if (!hz_) {
        for (int k = 0; k != nocc; ++k)
          mo1e(ij) -= 0.5*buf2e->element(k+i*nocc, k+j*nocc);
      }
    }
  }
}


void MOFile::update_1ext_ints(const shared_ptr<const Matrix>& coeff) {
  assert(mo2e_1ext_->nocc() == nocc_);
  assert(coeff->ndim() == nocc_);
  mo2e_1ext_ = mo2e_1ext_->transform_occ(coeff);
}


shared_ptr<const Matrix> Jop::compute_mo1e(const int nstart, const int nfence) {
  const MatView ocoeff = coeff_->slice(nstart, nfence);
  return make_shared<const Matrix>(ocoeff % *core_fock_ * ocoeff);
}


shared_ptr<const Matrix> Jop::compute_mo2e(const int nstart, const int nfence) {

  assert(nfence-nstart > 0);
  const MatView cdata = coeff_->slice(nstart, nfence);

  // first half transformation
  shared_ptr<DFHalfDist> half = geom_->df()->compute_half_transform(cdata);

  // second index transformation and (D|ii) = J^-1/2_DE (E|ii)
  // TODO : DFDistT needs to be modified to handle cases where number of nodes is larger than half->nocc() * cdata.mdim()
  shared_ptr<DFFullDist> buf;
  if (half->nocc() * cdata.mdim() > mpi__->size()) {
    buf = half->compute_second_transform(cdata)->apply_J();
  } else {
    buf = half->apply_J()->compute_second_transform(cdata);
  }

  // we want to store half-transformed quantity for latter convenience
  mo2e_1ext_ = half;

  // assembles (ii|ii) = (ii|D)(D|ii)
  return buf->form_4index(buf, 1.0);

}


