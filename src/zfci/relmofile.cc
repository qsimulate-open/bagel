//
// BAGEL - Parallel electron correlation program.
// Filename: relmofile.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <src/util/f77.h>
#include <src/zfci/relmofile.h>
#include <src/rel/dfock.h>
#include <src/rel/reldffull.h>

using namespace std;
using namespace bagel;

RelMOFile::RelMOFile(const shared_ptr<const Reference> ref, const string method)
: ZMOFile_Base(ref, method) {

  do_df_ = geom_->df().get();
  if (!do_df_) throw runtime_error("for the time being I gave up maintaining non-DF codes.");
}


double RelMOFile::create_Jiiii(const int nstart, const int nfence) {
  // first compute all the AO integrals in core
  norb_rel_ = nfence - nstart;
  nbasis_ = geom_->nbasis();
  geom_ = geom_->relativistic(false);
  assert(dynamic_pointer_cast<const RelReference>(ref_));
  auto relref = dynamic_pointer_cast<const RelReference>(ref_);

  // tempolary integral files
  shared_ptr<const ZMatrix> buf1e, buf2e;

  // one electron part
  double core_energy = 0;
  tie(buf1e, core_energy) = compute_mo1e(nstart, nfence);

  // two electron part.
  // this fills mo2e_1ext_ and returns buf2e which is an ii/ii quantity
  buf2e = compute_mo2e(nstart, nfence);

  // compute the unitary matrix that block-diagonalizes integrals 
  shared_ptr<const ZMatrix> umat = kramers();
  // set subblocks to mo1e_ and mo2e_
  tie(buf1e, buf2e) = kramers_block_diagonalize(umat, buf1e, buf2e);

  compress(buf1e, buf2e);
  return core_energy;
}


shared_ptr<const ZMatrix> RelMOFile::kramers() const {
  throw runtime_error("to be implemented");
  return shared_ptr<const ZMatrix>();
}


tuple<shared_ptr<const ZMatrix>,shared_ptr<const ZMatrix>>
 RelMOFile::kramers_block_diagonalize(shared_ptr<const ZMatrix> umat,shared_ptr<const ZMatrix> buf1e, shared_ptr<const ZMatrix> buf2e) {

  throw runtime_error("to be implemented");
}


//TODO input matrices are now in block diagonal form and the sizes must be checked
void RelMOFile::compress(shared_ptr<const ZMatrix> buf1e, shared_ptr<const ZMatrix> buf2e) {

  const int nocc = norb_rel_;
  sizeij_ = nocc*nocc;
  mo2e_ = unique_ptr<complex<double>[]>(new complex<double>[sizeij_*sizeij_]);
  copy_n(buf2e->data(), sizeij_*sizeij_, mo2e_.get());

  // h'ij = hij - 0.5 sum_k (ik|kj)
  const int size1e = nocc*nocc;
  mo1e_ = unique_ptr<complex<double>[]>(new complex<double>[size1e]);
  int ij = 0;
  for (int i = 0; i != nocc; ++i) {
    for (int j = 0; j != nocc; ++j, ++ij) {
      mo1e_[ij] = buf1e->element(j,i);
      for (int k = 0; k != nocc; ++k)
        mo1e_[ij] -= 0.5*buf2e->data(j+k*nocc+k*nocc*nocc+i*nocc*nocc*nocc);
    }
  }
}


tuple<shared_ptr<const ZMatrix>, double> RelJop::compute_mo1e(const int nstart, const int nfence) {

  auto relref = dynamic_pointer_cast<const RelReference>(ref_);
  complex<double> core_energy = 0.0;
  core_energy = 1e100;

  shared_ptr<ZMatrix> dfock0 = make_shared<RelHcore>(geom_);
  //TODO rhf density matrix as in ZMOFile, may not be implemented yet for dfock?
#if 0
  shared_ptr<RelHcore> relhcore = make_shared<RelHcore>(geom_);
  if (nstart != 0) {
    dfock0 = make_shared<DFock>(geom_, relhcore, relref->relcoeff(),true,true,true);
    core_energy = (*** RHF DENSITY MATRIX HERE *** (*relhcore+*dfock0)).trace() * 0.5;
  }
  dfock0->fill_upper();
#endif

  // Hij = relcoeff(T) * relhcore * relcoeff
  shared_ptr<const ZMatrix> coeff = relref->relcoeff()->slice(nstart, nfence);
  core_dfock_ = make_shared<ZMatrix>(*coeff % *dfock0 * *coeff);

  assert(fabs(core_energy.imag())<1e-10);
  return make_tuple(core_dfock_, core_energy.real());
}


shared_ptr<const ZMatrix> RelJop::compute_mo2e(const int nstart, const int nfence) {
//slightly modified code from rel/dmp2.cc to form 3 index integrals that we can build into 4 index with form4index
  const size_t norb_rel_ = nfence - nstart;
  if (norb_rel_ < 1) throw runtime_error("no correlated electrons");

  auto relref = dynamic_pointer_cast<const RelReference>(ref_);
  assert(geom_->nbasis()*4 == relref->relcoeff()->ndim());
  assert(geom_->nbasis()*2 == relref->relcoeff()->mdim());

  // Separate Coefficients into real and imaginary
  // correlated occupied orbitals
  array<shared_ptr<const Matrix>, 4> rocoeff;
  array<shared_ptr<const Matrix>, 4> iocoeff;
  for (int i = 0; i != 4; ++i) {
    const size_t nbasis = geom_->nbasis();
    shared_ptr<const ZMatrix> oc = relref->relcoeff()->get_submatrix(i*nbasis, nstart, nbasis, nfence);
    rocoeff[i] = oc->get_real_part();
    iocoeff[i] = oc->get_imag_part();
  }

  // (1) make DFDists
  vector<shared_ptr<const DFDist>> dfs;
  dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) first-transform
  list<shared_ptr<RelDFHalf>> half_complex = DFock::make_half_complex(dfdists, rocoeff, iocoeff);
  for (auto& i : half_complex)
    i = i->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_exch;
  for (auto& i : half_complex) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
  }
  half_complex.clear();
  DFock::factorize(half_complex_exch);

  // (4) compute (gamma|ii)
  list<shared_ptr<RelDFFull>> dffull;
  for (auto& i : half_complex_exch)
    dffull.push_back(make_shared<RelDFFull>(i, rocoeff, iocoeff));
  DFock::factorize(dffull);
  dffull.front()->scale(dffull.front()->fac()); // take care of the factor
  assert(dffull.size() == 1);
  shared_ptr<const RelDFFull> full = dffull.front();

  shared_ptr<const ZMatrix> buf = full->form_4index(full, 1.0);

  return buf;
}
