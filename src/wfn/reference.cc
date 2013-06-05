//
// BAGEL - Parallel electron correlation program.
// Filename: reference.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <boost/regex.hpp>

#include <src/util/f77.h>
#include <src/wfn/reference.h>
#include <src/fci/knowles.h>
#include <src/osint/overlapbatch.h>
#include <src/util/mixedbasis.h>
#include <src/util/lexical_cast.h>

using namespace std;
using namespace bagel;

Reference::Reference(shared_ptr<const Geometry> g, shared_ptr<const Coeff> c,
                     const int _nclosed, const int _nact, const int _nvirt,
                     const double en,
                     const vector<shared_ptr<RDM<1>>>& _rdm1, const vector<shared_ptr<RDM<2>>>& _rdm2,
                     shared_ptr<const RDM<1>> _rdm1_av, shared_ptr<const RDM<2>> _rdm2_av)
 : geom_(g), coeff_(c), energy_(en), hcore_(make_shared<Hcore>(geom_)), nclosed_(_nclosed), nact_(_nact), nvirt_(_nvirt), nstate_(1), rdm1_(_rdm1), rdm2_(_rdm2),
   rdm1_av_(_rdm1_av), rdm2_av_(_rdm2_av) {

  // we need to make sure that all the quantities are consistent in every MPI process
  mpi__->broadcast(coeff_->data(), coeff_->size(), 0);
  for (auto& i : rdm1_)
    mpi__->broadcast(i->data(), i->size(), 0);
  for (auto& i : rdm2_)
    mpi__->broadcast(i->data(), i->size(), 0);
  if (rdm1_av_)
    mpi__->broadcast_force(rdm1_av_->data(), rdm1_av_->size(), 0);
  if (rdm2_av_)
    mpi__->broadcast_force(rdm2_av_->data(), rdm2_av_->size(), 0);

  //if (nact_ && rdm1_.empty())
  //  throw logic_error("If nact != 0, Reference::Reference wants to have RDMs.");

}


shared_ptr<Matrix> Reference::rdm1_mat(shared_ptr<const RDM<1>> active) const {
  if (nact_)
    return active->rdm1_mat(geom_, nclosed_);
  else {
    auto out = make_shared<Matrix>(nocc(), nocc());
    for (int i = 0; i != nclosed_; ++i) out->element(i,i) = 2.0;
    return out;
  }
}


shared_ptr<Dvec> Reference::civectors() const {
  // Default to HarrisonZarrabian method
  shared_ptr<FCI> fci = make_shared<KnowlesHandy>(boost::property_tree::ptree(), shared_from_this(), nclosed_, nact_, nstate_);
  fci->compute();
  return fci->civectors();
}


// TODO this is a very bad implementation, since it recomputes FCI; should be replaced in somewhere.
tuple<shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> Reference::compute_rdm34(const int i) const {
  // Default to HarrisonZarrabian method
  shared_ptr<FCI> fci = make_shared<KnowlesHandy>(boost::property_tree::ptree(), shared_from_this(), nclosed_, nact_, nstate_);
  fci->compute();
  fci->compute_rdm12();
  return fci->compute_rdm34(i);
}


shared_ptr<const Reference> Reference::project_coeff(shared_ptr<const Geometry> geomin) const {
  shared_ptr<Matrix> snew = make_shared<Overlap>(geomin);
  snew->inverse_symmetric();
  MixedBasis<OverlapBatch> mixed(geom_, geomin);
  auto c = make_shared<Coeff>(*snew * mixed * *coeff_);

  return make_shared<Reference>(geomin, c, nclosed_, nact_, geomin->nbasis()-nclosed_-nact_, energy_);
}


void Reference::set_eig(const std::vector<double>& eig) {
  eig_ = eig;
  mpi__->broadcast(&eig_[0], eig_.size(), 0);
}


void Reference::set_erdm1(const shared_ptr<const Matrix> o) {
  mpi__->broadcast(o->data(), o->size(), 0);
  erdm1_ = o;
}


void Reference::set_coeff_AB(const shared_ptr<const Coeff> a, const shared_ptr<const Coeff> b) {
  mpi__->broadcast(a->data(), a->size(), 0);
  mpi__->broadcast(b->data(), b->size(), 0);
  coeffA_ = a;
  coeffB_ = b;
}

// This function currently assumes it is being called on a Reference object with no defined active space
shared_ptr<const Reference> Reference::set_active(set<int> active_indices) const {
  const int nbasis = geom_->nbasis();

  int nactive = active_indices.size();

  int nclosed = nclosed_;
  int nvirt = nbasis - nclosed;
  for (auto& iter : active_indices) {
    if (iter < nclosed_) --nclosed;
    else --nvirt;
  }

  auto tmp_coeff = make_shared<Matrix>(nbasis, nbasis);

  int iclosed = 0;
  int iactive = nclosed;
  int ivirt = nclosed + nactive;
  for (int i = 0; i < nclosed_; ++i) {
    if ( active_indices.find(i) == active_indices.end() ) {
      copy_n(coeff_->element_ptr(0,i), nbasis, tmp_coeff->element_ptr(0,iclosed));
      ++iclosed;
    }
    else {
      copy_n(coeff_->element_ptr(0,i), nbasis, tmp_coeff->element_ptr(0,iactive));
      ++iactive;
    }
  }

  for (int i = nclosed_; i < nbasis; ++i) {
    if ( active_indices.find(i) == active_indices.end() ) {
      copy_n(coeff_->element_ptr(0,i), nbasis, tmp_coeff->element_ptr(0,ivirt));
      ++ivirt;
    }
    else {
      copy_n(coeff_->element_ptr(0,i), nbasis, tmp_coeff->element_ptr(0,iactive));
      ++iactive;
    }
  }

  auto out_coeff = make_shared<const Coeff>(*tmp_coeff);
  return make_shared<Reference>(geom_, out_coeff, nclosed, nactive, nvirt);
}

shared_ptr<const Reference> Reference::set_active(string active_string) const {
  boost::regex r("(\\d+)");
  boost::smatch what;

  auto start = active_string.cbegin();
  auto end = active_string.cend();

  set<int> active_set;

  while( boost::regex_search(start, end, what, r) ) {
    string int_string(what[1].first, what[1].second);
    active_set.insert(lexical_cast<int>(int_string) - 1);
    start = what[0].second;
  }

  return set_active(active_set);
}
