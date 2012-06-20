//
// Newint - Parallel electron correlation program.
// Filename: supercigrad.cc
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


#include <src/grad/gradeval.h>
#include <src/grad/cpcasscf.h>
#include <src/casscf/supercigrad.h>
#include <src/util/pairfile.h>

using namespace std;

template<typename T>
static string tostring(const T i) {
  stringstream ss;
  ss << i; 
  return ss.str();
};


template<>
std::shared_ptr<GradFile> GradEval<SuperCIGrad>::compute() {

  // related to denominators
  const int nbasis = ref_->geom()->nbasis();
  // TODO TODO TODO temp... to have it compiled
  vector<double> eig(nbasis*nbasis);

  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();

  // TODO this is redundant, though...
  shared_ptr<DF_Half> half = ref_->geom()->df()->compute_half_transform(ref_->coeff()->data(), nocc);

  shared_ptr<FCI> fci;
  {
    multimap<string, string> fullci_data;
    fullci_data.insert(make_pair("ncore", tostring<int>(nclosed)));
    fullci_data.insert(make_pair("norb",  tostring<int>(nact)));
    fci = shared_ptr<FCI>(new FCI(fullci_data, ref_->geom(), ref_));
  }

  int la, lb; tie(la, lb) = fci->len_string();
  shared_ptr<Matrix1e> g0(new Matrix1e(ref_->geom())); 
  shared_ptr<Dvec> g1(new Dvec(lb, la, ref_->nstate()));
  shared_ptr<PairFile<Matrix1e, Dvec> > grad(new PairFile<Matrix1e, Dvec>(g0, g1));
  shared_ptr<CPCASSCF> cp(new CPCASSCF(grad, eig, half, ref_, fci));

  shared_ptr<PairFile<Matrix1e, Dvec> > zvec = cp->solve();

  std::shared_ptr<GradFile> out(new GradFile(3*geom_->natom()));
  return out;
}

