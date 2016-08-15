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

void Box::init() {

  extent_ = 0;
  if (rank_ == 0) {
    for (auto& i : sp_) {
      centre_[0] += i->centre(0);
      centre_[1] += i->centre(1);
      centre_[2] += i->centre(2);
      if (extent_ < i->extent()) extent_ = i->extent();
    }
    centre_[0] /= nsp();
    centre_[1] /= nsp();
    centre_[2] /= nsp();
  } else {
    assert (!child_.empty());
    for (int i = 0; i != nchild(); ++i) {
      shared_ptr<const Box> c = child(i);
      centre_[0] += c->centre(0);
      centre_[1] += c->centre(1);
      centre_[2] += c->centre(2);
      if (extent_ < c->extent()) extent_ = c->extent();
    }
    centre_[0] /= nchild();
    centre_[1] /= nchild();
    centre_[2] /= nchild();
  }
}


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


void Box::get_neigh(vector<shared_ptr<Box>> box, const int ws) {

  neigh_.resize(box.size());
  int nn = 0;
  for (auto& b : box) {
    if (b->rank() == rank_) {
      if (is_neigh(b)) {
        neigh_[nn] = b;
        ++nn;
      }
    }
  }
  neigh_.resize(nn);
}


bool Box::is_neigh(shared_ptr<const Box> box, const int ws) const {

  double rr = 0;
  for (int i = 0; i != 3; ++i)
    rr += pow(centre_[i] - box->centre(i), 2);

  const bool out = (sqrt(rr) > (1+ws)*(extent_ + box->extent()));
  return out;
}
