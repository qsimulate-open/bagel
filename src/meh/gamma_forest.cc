//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_forest.cc
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <utility>

#include <src/meh/gamma_forest.h>

using namespace bagel;
using namespace std;

GammaTree::GammaTree(shared_ptr<const Dvec> ket) : ket_(ket) {
  base_ = make_shared<GammaBranch>();
  const int nops = 4;

  for (int i = 0; i < nops; ++i) {
    base_->branch(i) = make_shared<GammaBranch>();
    for (int j = 0; j < nops; ++j) {
      base_->branch(i)->branch(j) = make_shared<GammaBranch>();
      for (int k = 0; k < nops; ++k) {
        base_->branch(i)->branch(j)->branch(k) = make_shared<GammaBranch>();
      }
    }
  }
}
