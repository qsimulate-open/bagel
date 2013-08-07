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

// TODO one line change in order to add a function to block diagonalize
double RelMOFile::create_Jiiii(const int nstart, const int nfence) {
  // first compute all the AO integrals in core
  nocc_ = nfence - nstart;
  nbasis_ = geom_->nbasis();
  const int nbasis = nbasis_;
  relgeom_ = geom_->relativistic(false);
  relref = dynamic_pointer_cast<const RelReference>(ref_);

  // one electron part
  double core_energy = 0;
  shared_ptr<const ZMatrix> buf1e;
  tie(buf1e, core_energy) = compute_mo1e(nstart, nfence);

  // two electron part.
  // this fills mo2e_1ext_ and returns buf2e which is an ii/ii quantity
  unique_ptr<complex<double>[]> buf2e = compute_mo2e(nstart, nfence);

  // TODO block diaganolize according kramers symmetry
  //kramers_block(buf1e, buf2e);

  compress(buf1e, buf2e);
  return core_energy;
}

void RelMOFile::kramers_block(shared_ptr<const ZMatrix> buf1e, unique_ptr<complex<double>[]>& buf2e) {
//
//block diagonalize 1e and 2e matrices according to kramers symmetry by a unitary operator
//
//1) sort degenerate orbitals into a vector
//2) loop over vector and apply projection
//
}

//does not need to be changed as long as kramers_block is in the same format
void RelMOFile::compress(shared_ptr<const ZMatrix> buf1e, unique_ptr<complex<double>[]>& buf2e) {

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

//
//Relativistic 1 and 2e functions in progress
//
tuple<shared_ptr<const ZMatrix>, double> RelJop::compute_mo1e(const int nstart, const int nfence) {

  const size_t nbasis = geom_->nbasis();

  complex<double> core_energy = 0.0;

  auto relhcore = make_shared<RelHcore>(relgeom_);

  shared_ptr<const ZMatrix> coeff = relref->relcoeff()->slice(nstart, nfence);

  // Hij = relcoeff(T) * relhcore * relcoeff
  core_dfock_ = make_shared<ZMatrix>(*coeff % *relhcore * *coeff); 

  //TODO include some density adjustment? see zmofile
//core_energy = relhcore->trace();
  // FIXME
  core_energy = 1.0e100;
  assert(fabs(core_energy.imag())<1e-10);

  return make_tuple(core_dfock_, core_energy.real());
}


unique_ptr<complex<double>[]> RelJop::compute_mo2e(const int nstart, const int nfence) {
//slightly modified code from rel/dmp2.cc to form 3 index integrals that we can build into 4 index with form4index
  const size_t nocc = nfence - nstart; 
  if (nocc < 1) throw runtime_error("no correlated electrons");

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
  dfs = relgeom_->dfs()->split_blocks();
  dfs.push_back(relgeom_->df());
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
  cout << "    * 3-index integral transformation done" << endl;

  unique_ptr<complex<double>[]> out(new complex<double>[nocc*nocc*nocc*nocc]);
  {
    shared_ptr<const ZMatrix> buf = full->form_4index(full, 1.0);
    copy_n(buf->data(), nocc*nocc*nocc*nocc, out.get());
  }
  return out;
}
