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
#include <src/fci/fci.h>
#include <src/osint/overlapbatch.h>
#include <src/util/mixedbasis.h>

using namespace std;

Reference::Reference(shared_ptr<const Geometry> g, shared_ptr<const Coeff> c,
                     const int _nclosed, const int _nact, const int _nvirt,
                     const double en,
                     const vector<shared_ptr<RDM<1> > >& _rdm1, const vector<shared_ptr<RDM<2> > >& _rdm2,
                     shared_ptr<const RDM<1> > _rdm1_av, shared_ptr<const RDM<2> > _rdm2_av)
 : geom_(g), coeff_(c), energy_(en), hcore_(new Hcore(geom_)), nclosed_(_nclosed), nact_(_nact), nvirt_(_nvirt), nstate_(1), rdm1_(_rdm1), rdm2_(_rdm2),
   rdm1_av_(_rdm1_av), rdm2_av_(_rdm2_av) {

  if (nact_ && rdm1_.empty())
    throw logic_error("If nact != 0, Reference::Reference wants to have RDMs.");

}


shared_ptr<Matrix1e> Reference::rdm1_mat(shared_ptr<const RDM<1> > active) const {
  if (nact_)
    return active->rdm1_mat(geom_, nclosed_);
  else {
    shared_ptr<Matrix1e> out(new Matrix1e(geom_, nocc(), nocc()));
    for (int i = 0; i != nclosed_; ++i) out->element(i,i) = 2.0;
    return out;
  }
}


shared_ptr<Dvec> Reference::civectors() const {
  shared_ptr<FCI> fci(new FCI(multimap<string, string>(), shared_from_this(), nclosed_, nact_, nstate_));
  fci->compute();
  return fci->civectors();
}


shared_ptr<const Coeff> Reference::project_coeff(shared_ptr<const Geometry> geomin) const {
  shared_ptr<const Coeff> out;

  if (*geom_ == *geomin) {
    out = coeff_;
  } else {
    // in this case we first form overlap matrices
    shared_ptr<Overlap> snew(new Overlap(geomin));
    shared_ptr<Overlap> sold(new Overlap(geom_));
    shared_ptr<MixedBasis<OverlapBatch> > mixed(new MixedBasis<OverlapBatch>(geomin, geom_));
mixed->print("mixed", 24);
  }

  return out; 
}
