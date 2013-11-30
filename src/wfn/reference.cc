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

#include <src/wfn/reference.h>
#include <src/fci/knowles.h>
#include <src/integral/os/overlapbatch.h>
#include <src/molecule/mixedbasis.h>

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
  if (coeff_)
    mpi__->broadcast(const_pointer_cast<Coeff>(coeff_)->data(), coeff_->size(), 0);

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
    return active->rdm1_mat(nclosed_);
  else {
    auto out = make_shared<Matrix>(nocc(), nocc());
    for (int i = 0; i != nclosed_; ++i) out->element(i,i) = 2.0;
    return out;
  }
}


// TODO should be updated to remove redundant FCI iterations
shared_ptr<Dvec> Reference::civectors() const {
  // Default to HarrisonZarrabian method
  shared_ptr<FCI> fci = make_shared<KnowlesHandy>(make_shared<const PTree>(), geom_, shared_from_this(), nclosed_, nact_, nstate_);
  fci->compute();
  return fci->civectors();
}


// TODO should be updated to remove redundant FCI iterations
shared_ptr<Dvec> Reference::rdm1deriv() const {
  shared_ptr<FCI> fci = make_shared<KnowlesHandy>(make_shared<const PTree>(), geom_, shared_from_this(), nclosed_, nact_, nstate_);
  fci->compute();
  shared_ptr<Dvec> out = fci->rdm1deriv();
  return out;
}


// TODO should be updated to remove redundant FCI iterations
shared_ptr<Dvec> Reference::rdm2deriv() const {
  shared_ptr<FCI> fci = make_shared<KnowlesHandy>(make_shared<const PTree>(), geom_, shared_from_this(), nclosed_, nact_, nstate_);
  fci->compute();
  shared_ptr<Dvec> out = fci->rdm2deriv();
  return out;
}


// TODO should be updated to remove redundant FCI iterations
shared_ptr<Dvec> Reference::rdm3deriv() const {
  shared_ptr<FCI> fci = make_shared<KnowlesHandy>(make_shared<const PTree>(), geom_, shared_from_this(), nclosed_, nact_, nstate_);
  fci->compute();
  shared_ptr<Dvec> out = fci->rdm3deriv();
  return out;
}


// TODO should be updated to remove redundant FCI iterations
shared_ptr<Dvec> Reference::rdm4deriv() const {
  shared_ptr<FCI> fci = make_shared<KnowlesHandy>(make_shared<const PTree>(), geom_, shared_from_this(), nclosed_, nact_, nstate_);
  fci->compute();
  shared_ptr<Dvec> out = fci->rdm4deriv();
  return out;
}


// TODO this is a very bad implementation, since it recomputes FCI; should be replaced in somewhere.
tuple<shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> Reference::compute_rdm34(const int i) const {
  // Default to HarrisonZarrabian method
  shared_ptr<FCI> fci = make_shared<KnowlesHandy>(make_shared<const PTree>(), geom_, shared_from_this(), nclosed_, nact_, nstate_);
  fci->compute();
  fci->compute_rdm12();
  return fci->compute_rdm34(i);
}


shared_ptr<const Reference> Reference::project_coeff(shared_ptr<const Geometry> geomin) const {
  shared_ptr<Matrix> snew = make_shared<Overlap>(geomin);
  snew->inverse_symmetric();
  MixedBasis<OverlapBatch> mixed(geom_, geomin);
  auto c = make_shared<Coeff>(*snew * mixed * *coeff_);

  auto out = make_shared<Reference>(geomin, c, nclosed_, nact_, geomin->nbasis()-nclosed_-nact_, energy_);
  if (coeffA_) {
    assert(coeffB_);
    out->coeffA_ = make_shared<Coeff>(*snew * mixed * *coeffA_);
    out->coeffB_ = make_shared<Coeff>(*snew * mixed * *coeffB_);
  }
  return out;
}


void Reference::set_eig(const std::vector<double>& eig) {
  eig_ = eig;
  mpi__->broadcast(&eig_[0], eig_.size(), 0);
}


void Reference::set_erdm1(const shared_ptr<const Matrix> o) {
  mpi__->broadcast(const_pointer_cast<Matrix>(o)->data(), o->size(), 0);
  erdm1_ = o;
}


void Reference::set_coeff_AB(const shared_ptr<const Coeff> a, const shared_ptr<const Coeff> b) {
  mpi__->broadcast(const_pointer_cast<Coeff>(a)->data(), a->size(), 0);
  mpi__->broadcast(const_pointer_cast<Coeff>(b)->data(), b->size(), 0);
  coeffA_ = a;
  coeffB_ = b;
}

// This function currently assumes it is being called on a Reference object with no defined active space
shared_ptr<Reference> Reference::set_active(set<int> active_indices) const {
  if (!coeff_) throw logic_error("Reference::set_active is not implemented for relativistic cases");
  const int naobasis = geom_->nbasis();
  const int nmobasis = coeff_->mdim();

  int nactive = active_indices.size();

  int nclosed = nclosed_;
  int nvirt = nmobasis - nclosed;
  for (auto& iter : active_indices) {
    if (iter < nclosed_) --nclosed;
    else --nvirt;
  }

  auto coeff = coeff_;
  auto tmp_coeff = make_shared<Matrix>(naobasis, nmobasis);

  int iclosed = 0;
  int iactive = nclosed;
  int ivirt = nclosed + nactive;

  auto cp = [&tmp_coeff, &naobasis, &coeff] (const int i, int& pos) { copy_n(coeff->element_ptr(0,i), naobasis, tmp_coeff->element_ptr(0, pos)); ++pos; };

  for (int i = 0; i < nmobasis; ++i) {
    if ( active_indices.find(i) != active_indices.end() ) cp(i, iactive);
    else if ( i < nclosed_ ) cp(i, iclosed);
    else cp(i, ivirt);
  }

  return make_shared<Reference>(geom_, make_shared<const Coeff>(*tmp_coeff), nclosed, nactive, nvirt);
}

// This function currently assumes it is being called on a Reference object with no defined active space
shared_ptr<Reference> Reference::set_ractive(set<int> ras1, set<int> ras2, set<int> ras3) const {
  if (!coeff_) throw logic_error("Reference::set_active is not implemented for relativistic cases");
  const int naobasis = geom_->nbasis();
  const int nmobasis = coeff_->mdim();

  const int nras1 = ras1.size();
  const int nras2 = ras2.size();
  const int nras3 = ras3.size();

  int nactive = nras1 + nras2 + nras3;

  set<int> total_active;
  total_active.insert(ras1.begin(), ras1.end());
  total_active.insert(ras2.begin(), ras2.end());
  total_active.insert(ras3.begin(), ras3.end());
  if (total_active.size() != nactive) throw runtime_error("Each orbital can occur in only one of RAS1, RAS2, or RAS3.");

  int nclosed = nclosed_;
  int nvirt = nmobasis - nclosed;
  for (auto& iter : total_active) {
    if (iter < nclosed_) --nclosed;
    else --nvirt;
  }

  auto coeff = coeff_;
  auto tmp_coeff = make_shared<Matrix>(naobasis, nmobasis);

  int iclosed = 0;
  int iras1 = nclosed;
  int iras2 = iras1 + nras1;
  int iras3 = iras2 + nras2;
  int ivirt = nclosed + nactive;

  auto cp = [&tmp_coeff, &naobasis, &coeff] (const int i, int& pos) { copy_n(coeff->element_ptr(0,i), naobasis, tmp_coeff->element_ptr(0, pos)); ++pos; };

  for (int i = 0; i < nmobasis; ++i) {
    if ( total_active.find(i) != total_active.end() ) {
      if (ras1.find(i) != ras1.end()) cp(i, iras1);
      else if (ras2.find(i) != ras2.end()) cp(i, iras2);
      else if (ras3.find(i) != ras3.end()) cp(i, iras3);
      else assert(false);
    }
    else if (i < nclosed_) cp(i, iclosed);
    else cp(i, ivirt);
  }

  return make_shared<Reference>(geom_, make_shared<const Coeff>(*tmp_coeff), nclosed, nactive, nvirt);
}
