//
// BAGEL - Parallel electron correlation program.
// Filename: breit2index.cc
// Copyright (C) 2013 Matthew Kelley
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
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


#include <stddef.h>
#include <src/rel/breit2index.h>

using namespace std;
using namespace bagel;

Breit2Index::Breit2Index(pair<const int, const int> index, shared_ptr<const Matrix> breit, shared_ptr<const Matrix> dat2)
 : index_(index), data_(make_shared<Matrix>(*dat2 % *breit * *dat2)) {
}


shared_ptr<Breit2Index> Breit2Index::cross() const {
  int i = index_.first;
  int j = index_.second;
  return make_shared<Breit2Index>(make_pair(j,i), data_);
}


void Breit2Index::print() const {
  cout << " Breit2Index::data_" << endl;
  data_->print();
}
