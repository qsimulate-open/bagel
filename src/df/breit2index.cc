//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: breit2index.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
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


#include <stddef.h>
#include <src/df/breit2index.h>

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
