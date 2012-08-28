//
// BAGEL - Parallel electron correlation program.
// Filename: fci.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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
#include <stdexcept>
#include <src/fci/fci.h>

using namespace std;
using namespace bagel;

FCI::FCI(std::multimap<std::string, std::string> idat, shared_ptr<const Reference> r, const int ncore, const int norb, const int nstate)
 : idata_(idat), ref_(r), geom_(r->geom()), ncore_(ncore), norb_(norb), nstate_(nstate) {
  common_init();
}


void FCI::common_init() {
  print_header();

  const bool frozen = read_input<bool>(idata_, "frozen", false);
  max_iter_ = read_input<int>(idata_, "maxiter", 100);
  max_iter_ = read_input<int>(idata_, "maxiter_fci", max_iter_);
  thresh_ = read_input<double>(idata_, "thresh", 1.0e-20);
  thresh_ = read_input<double>(idata_, "thresh_fci", thresh_);

  if (nstate_ < 0) nstate_ = read_input<int>(idata_, "nstate", 1);
  if (ncore_ < 0) ncore_ = read_input<int>(idata_, "ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (norb_  < 0) norb_ = read_input<int>(idata_, "norb", ref_->coeff()->ndim()-ncore_);

  // nspin is #unpaired electron 0:singlet, 1:doublet, 2:triplet, ... (i.e., Molpro convention).
  const int nspin = read_input<int>(idata_, "nspin", 0);
  if ((geom_->nele()+nspin) % 2 != 0) throw runtime_error("Invalid nspin specified");
  nelea_ = (geom_->nele()+nspin)/2 - ncore_;
  neleb_ = (geom_->nele()-nspin)/2 - ncore_;

  // TODO allow for zero electron (quick return)
  if (nelea_ <= 0 || neleb_ <= 0) throw runtime_error("#electrons cannot be zero/negative in FCI");
  for (int i = 0; i != nstate_; ++i) weight_.push_back(1.0/static_cast<double>(nstate_));

  // resizing rdm vectors (with null pointers)
  rdm1_.resize(nstate_);
  rdm2_.resize(nstate_);
  energy_.resize(nstate_);

  // construct a determinant space in which this FCI will be performed.
  det_ = shared_ptr<const Determinants>(new Determinants(norb_, nelea_, neleb_));

  // forms MO integrals and denominators.
  update(ref_->coeff());

}

FCI::~FCI() {

}



void FCI::print_timing_(const string label, int& time, std::vector<pair<string, double> >& timing) const {
  timing.push_back(make_pair(label, (::clock()-time)/static_cast<double>(CLOCKS_PER_SEC)));
  time = ::clock();
}
