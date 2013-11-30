//
// BAGEL - Parallel electron correlation program.
// Filename: mofile.cc
// Copyright (C) 2011 Toru Shiozaki
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

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <src/util/f77.h>
#include <src/fci/mofile.h>
#include <src/scf/scf.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;


MOFile::MOFile(const shared_ptr<const Reference> ref, const string method)
: geom_(ref->geom()), ref_(ref), coeff_(ref_->coeff()) {

  do_df_ = geom_->df().get();
  if (!do_df_) throw runtime_error("for the time being I gave up maintaining non-DF codes.");

  hz_ = (method=="HZ");
}


MOFile::MOFile(const shared_ptr<const Reference> ref, const shared_ptr<const Coeff> c, const string method)
: hz_(false), geom_(ref->geom()), ref_(ref), coeff_(c) {

  do_df_ = geom_->df().get();
  if (!do_df_) throw runtime_error("for the time being I gave up maintaining non-DF codes.");

  hz_ = (method=="HZ");
}


void MOFile::init(const int nstart, const int nfence) {

  // first compute all the AO integrals in core
  nocc_ = nfence - nstart;

  // core energy is set here
  if (nstart != 0) {
    shared_ptr<const Matrix> den = coeff_->form_density_rhf(nstart);
    core_fock_ = make_shared<Fock<1>>(geom_, ref_->hcore(), den, ref_->schwarz());
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
        const int ijo = (j + i*nocc)*nocc*nocc;
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
    SMITH::sort_indices<0,2,1,3,0,1,1,1>(buf2e->data(), mo2e_->data(), nocc_, nocc_, nocc_, nocc_);
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
  mo2e_1ext_->rotate_occ(coeff);
}


shared_ptr<const Matrix> Jop::compute_mo1e(const int nstart, const int nfence) {
  shared_ptr<const Matrix> ocoeff = coeff_->slice(nstart, nfence);
  return make_shared<const Matrix>(*ocoeff % *core_fock_ * *ocoeff);
}


shared_ptr<const Matrix> Jop::compute_mo2e(const int nstart, const int nfence) {

  const int nocc = nfence - nstart;
  assert(nocc > 0);
  shared_ptr<const Matrix> cdata = coeff_->slice(nstart, nfence);

  // first half transformation
  shared_ptr<DFHalfDist> half = geom_->df()->compute_half_transform(cdata);

  // second index transformation and (D|ii) = J^-1/2_DE (E|ii)
  shared_ptr<DFFullDist> buf = half->compute_second_transform(cdata)->apply_J();

  // we want to store half-transformed quantity for latter convenience
  mo2e_1ext_size_ = nocc*geom_->df()->naux()*geom_->nbasis();
  mo2e_1ext_ = half;

  // assembles (ii|ii) = (ii|D)(D|ii)
  return buf->form_4index(buf, 1.0);

}


