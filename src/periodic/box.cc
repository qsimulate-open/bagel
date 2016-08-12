//
// BAGEL - Parallel electron correlation program.
// Filename: box.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/util/f77.h>
#include <src/periodic/box.h>
#include <src/periodic/multipolebatch.h>
#include <src/integral/os/overlapbatch.h>
#include <src/periodic/localexpansion.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;

void Box::insert_sp(vector<shared_ptr<const ShellPair>> sp) {

  const int nsp = sp_.size();
  sp_.resize(nsp + sp.size());
  for (int i = 0; i != sp.size(); ++i)
    sp_[i + nsp] = sp[i];

}


void Box::insert_child(shared_ptr<const Box> child) {

  const int nchild = child_.size();
  if (child) {
    child_.resize(nchild + 1);
    child_[nchild] = child;
  }
}
