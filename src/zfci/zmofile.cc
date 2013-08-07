//
// BAGEL - Parallel electron correlation program.
// Filename: zmofile.cc
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell  <caldwell@u.northwestern.edu>
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
#include <src/zfci/zmofile.h>
#include <src/math/algo.h>
#include <src/scf/scf.h>

using namespace std;
using namespace bagel;

ZMOFile::ZMOFile(const shared_ptr<const Reference> ref, const string method)
: ZMOFile_Base(ref, method), core_fock_(make_shared<Matrix>(geom_->nbasis(), geom_->nbasis())) {

  do_df_ = geom_->df().get();
  if (!do_df_) throw runtime_error("for the time being I gave up maintaining non-DF codes.");

  hz_ = (method=="HZ");
}


ZMOFile::ZMOFile(const shared_ptr<const Reference> ref, const shared_ptr<const Coeff> c, const string method)
: ZMOFile_Base(ref, method), core_fock_(make_shared<Matrix>(geom_->nbasis(), geom_->nbasis())), coeff_(c) {

  do_df_ = geom_->df().get();
  if (!do_df_) throw runtime_error("for the time being I gave up maintaining non-DF codes.");

  hz_ = (method=="HZ");
}


double ZMOFile::create_Jiiii(const int nstart, const int nfence) {
  // first compute all the AO integrals in core
  nocc_ = nfence - nstart;
  nbasis_ = geom_->nbasis();
  const int nbasis = nbasis_;

  // one electron part
  double core_energy;
  shared_ptr<const ZMatrix> buf1e;
  tie(buf1e, core_energy) = compute_mo1e(nstart, nfence);

  // two electron part.
  // this fills mo2e_1ext_ and returns buf2e which is an ii/ii quantity
  unique_ptr<complex<double>[]> buf2e = compute_mo2e(nstart, nfence);
  compress(buf1e, buf2e);
  return core_energy;
}

void ZMOFile::compress(shared_ptr<const ZMatrix> buf1e, unique_ptr<complex<double>[]>& buf2e) {

  const int nocc = nocc_;
  sizeij_ = nocc*nocc;
  mo2e_ = unique_ptr<complex<double>[]>(new complex<double>[sizeij_*sizeij_]);
  copy_n(buf2e.get(), sizeij_*sizeij_, mo2e_.get());

  // h'ij = hij - 0.5 sum_k (ik|kj)
  const int size1e = nocc*nocc;
  mo1e_ = unique_ptr<complex<double>[]>(new complex<double>[size1e]);
  int ij = 0;
  for (int i = 0; i != nocc; ++i) {
    for (int j = 0; j != nocc; ++j, ++ij) {
      mo1e_[ij] = buf1e->element(j,i);
      for (int k = 0; k != nocc; ++k)
        mo1e_[ij] -= 0.5*buf2e[j+k*nocc+k*nocc*nocc+i*nocc*nocc*nocc];
    }
  }
}
#if 0
void ZMOFile::update_1ext_ints(const shared_ptr<const Matrix>& coeff) {
  assert(mo2e_1ext_->nocc() == nocc_);
  assert(coeff->ndim() == nocc_);
  mo2e_1ext_->rotate_occ(coeff->data());
}
#endif

//NOTE because complex scf is not yet implemented this returns all zero imaginary terms
tuple<shared_ptr<const ZMatrix>, double> ZJop::compute_mo1e(const int nstart, const int nfence) {

  const int ncore = nstart;
  double core_energy = 0.0;

  auto fock0 = make_shared<Matrix>(*ref_->hcore());
  // if core fock operator is not the same as hcore...
  if (nstart != 0) {
    shared_ptr<Matrix> den = coeff_->form_density_rhf(ncore);
    fock0 = make_shared<Fock<1>>(geom_, ref_->hcore(), den, ref_->schwarz());
    core_energy = (*den * (*ref_->hcore()+*fock0)).trace() * 0.5;
    *core_fock_ = *fock0;
  }
  fock0->fill_upper();

  shared_ptr<const Matrix> ocoeff = coeff_->slice(nstart, nfence);
  auto mat = make_shared<const ZMatrix>(*ocoeff % *fock0 * *ocoeff, 1.0);

  return make_tuple(mat, core_energy);
}
//NOTE because complex df is not yet implemented this returns all zero imaginary terms
unique_ptr<complex<double>[]> ZJop::compute_mo2e(const int nstart, const int nfence) {

  const int nocc = nfence - nstart;
  assert(nocc > 0);
  shared_ptr<const Matrix> cdata = coeff_->slice(nstart, nfence);

  // first half transformation
  shared_ptr<DFHalfDist> half = geom_->df()->compute_half_transform(cdata);

  // second index transformation and (D|ii) = J^-1/2_DE (E|ii)
  shared_ptr<DFFullDist> buf = half->compute_second_transform(cdata)->apply_J();

  // we want to store half-transformed quantity for latter convenience
  mo2e_1ext_size_ = nocc*geom_->df()->naux()*nbasis_;
  mo2e_1ext_ = half;

  // assembles (ii|ii) = (ii|D)(D|ii)
  shared_ptr<const Matrix> med = buf->form_4index(buf, 1.0);
  unique_ptr<complex<double>[]> out(new complex<double>[nocc*nocc*nocc*nocc]);
  for (int i = 0; i != nocc*nocc*nocc*nocc; i++) out[i] = med->data(i);
  return out;
}

