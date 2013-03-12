//
// BAGEL - Parallel electron correlation program.
// Filename: dimer_levelshift.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#include <array>
#include <vector>
#include <set>
#include <string>
#include <memory>

#include <src/dimer/dimer_levelshift.h>
#include <src/scf/overlap.h>
#include <src/scf/coeff.h>

using namespace bagel;
using namespace std;

ShiftDimer::ShiftDimer(shared_ptr<const Dimer> dimer, const double shift_parameter) : LevelShift(dimer->nbasis().first, shift_parameter) {
  nbasis_ = dimer->nbasis();
  dimerbasis_ = dimer->dimerbasis();

  const int nstart = nbasis_.first;
  const int nfence = nstart + nbasis_.second;
  subspace_ = dimer->proj_coeff()->slice(nstart, nfence);

  S_ = shared_ptr<const Matrix>(new const Overlap(dimer->sgeom()));
  subspace_projector_ = shared_ptr<const Matrix>(new const Matrix( (*S_) * (*dimer->proj_coeff()) ));
}

void ShiftDimer::shift(Matrix& fock_mo, shared_ptr<const Coeff> coeff) {
  // Project onto subspace
  Matrix overlap( (*coeff) % (*subspace_projector_) );

  // multimap will order according to overlap
  multimap<double, int> Anorms, Bnorms;

  double* odata = overlap.data();

  const int nbasisA = nbasis_.first;
  const int nbasisB = nbasis_.second;

  for(int i = 0; i < dimerbasis_; ++i, ++odata) {
    double A_norm = ddot_(nbasisA, odata, dimerbasis_, odata, dimerbasis_);
    Anorms.insert(make_pair(A_norm, i));
  }

  odata = overlap.data() + dimerbasis_ * nbasisA;
  for(int i = 0; i < dimerbasis_; ++i, ++odata) {
    double B_norm = ddot_(nbasisB, odata, dimerbasis_, odata, dimerbasis_);
    Bnorms.insert(make_pair(B_norm, i));
  }

  set<int> Aset, Bset;

  auto aiter = Anorms.begin();
  auto biter = Bnorms.rbegin();
  for(int i = 0; i < nbasisB; ++i, ++aiter, ++biter) {
    Aset.insert(aiter->second);
    Bset.insert(biter->second);
  }

  if( Aset != Bset ) {
    cout << "Uhoh, we may be in trouble." << endl;
  }
  else {
    for (auto& aiter : Aset) {
      // Shift the highest
      fock_mo.element(aiter,aiter) += shift_parameter_;
    }
  }
}

void ShiftDimer::print_mo_data(shared_ptr<const Coeff> coeff) {
  // Project onto subspace
  Matrix overlap( (*coeff) % (*subspace_projector_) );

  // multimap will order according to overlap
  vector<double> Anorms, Bnorms;

  double* odata = overlap.data();

  const int nbasisA = nbasis_.first;
  const int nbasisB = nbasis_.second;

  for(int i = 0; i < dimerbasis_; ++i, ++odata) {
    double A_norm = ddot_(nbasisA, odata, dimerbasis_, odata, dimerbasis_);
    Anorms.push_back(A_norm);
  }

  odata = overlap.data() + dimerbasis_ * nbasisA;
  for(int i = 0; i < dimerbasis_; ++i, ++odata) {
    double B_norm = ddot_(nbasisB, odata, dimerbasis_, odata, dimerbasis_);
    Bnorms.push_back(B_norm);
  }

  cout << endl << "Subspace overlaps" << endl;
  cout << setw(5) << "MO" << setw(10) << "A" << setw(10) << "B" << endl;
  for(int i = 0; i < dimerbasis_; ++i)
    cout << setw(5) << i << setw(10) << setprecision(6) << Anorms.at(i) << setw(10) << Bnorms.at(i) << endl;
}
