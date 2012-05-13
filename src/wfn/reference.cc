//
// Newint - Parallel electron correlation program.
// Filename: reference.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/wfn/reference.h>

using namespace std;

Reference::Reference(shared_ptr<const Geometry> g, shared_ptr<Coeff> c,  const double en, shared_ptr<Hcore> h, const vector<double>& s,
                     const int& ncl, const int& nac, const int& nvi,
                     const vector<shared_ptr<RDM<1> > > _rdm1, const vector<shared_ptr<RDM<2> > > _rdm2)
 : geom_(g), coeff_(c), energy_(en), hcore_(h), schwarz_(s), nclosed_(ncl), nact_(nac), nvirt_(nvi), rdm1_(_rdm1), rdm2_(_rdm2) {

  if (nact_ && (rdm1_.empty() || rdm2_.empty()))
    throw logic_error("If nact != 0, Reference::Reference wants to have RDMs.");

}


shared_ptr<Matrix1e> Reference::rdm1() const {
  shared_ptr<Matrix1e> out(new Matrix1e(geom_, nocc(), nocc()));

  // first fill in diagonal elements for closed orbitals
  for (int i = 0; i != nclosed_; ++i) {
    out->element(i,i) = 2.0;
  }

  // TODO
  for (int i = 0; i != nact_; ++i) {
    throw logic_error("not yet implemented!!!!");
  }
  return out;
}
