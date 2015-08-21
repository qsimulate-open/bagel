//
// BAGEL - Parallel electron correlation program.
// Filename: reldffullt.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/df/reldffullt.h>

using namespace std;
using namespace bagel;

RelDFFullT::RelDFFullT(shared_ptr<const RelDFFull> full, shared_ptr<const StaticDist> dist) : basis_(full->basis().front()){
  assert(full->basis().size() == 1);
  dffull_[0] = make_shared<DFDistT>(full->get_real(), dist);
  dffull_[1] = make_shared<DFDistT>(full->get_imag(), dist);
}


void RelDFFullT::discard_df() {
  dffull_[0]->discard_df();
  dffull_[1]->discard_df();
}


vector<pair<shared_ptr<Matrix>,shared_ptr<Matrix>>> RelDFFullT::get_slice(const int start, const int end) const {
  assert(dffull_[0]->nblocks() == 1 && dffull_[1]->nblocks() == 1);
  return {make_pair(dffull_[0]->get_slice(start, end).at(0), dffull_[1]->get_slice(start, end).at(0))};
}


shared_ptr<ZMatrix> RelDFFullT::replicate() const {
  shared_ptr<const Matrix> r = dffull_[0]->replicate();
  shared_ptr<const Matrix> i = dffull_[1]->replicate();
  return make_shared<ZMatrix>(*r, *i);
}
