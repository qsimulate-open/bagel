//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reference.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/wfn/reference.h>
#include <src/wfn/relreference.h>
#include <src/integral/os/overlapbatch.h>
#include <src/mat1e/mixedbasis.h>
#include <src/ci/fci/fci.h>
#include <src/util/io/moldenin.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Reference)

using namespace std;
using namespace bagel;

Reference::Reference(shared_ptr<const Geometry> g, shared_ptr<const Coeff> c,
                     const int _nclosed, const int _nact, const int _nvirt,
                     const vector<double> en,
                     shared_ptr<const VecRDM<1>> _rdm1, shared_ptr<const VecRDM<2>> _rdm2,
                     shared_ptr<const RDM<1>> _rdm1_av, shared_ptr<const RDM<2>> _rdm2_av,
                     shared_ptr<const CIWfn> ci)
 : geom_(g), noccA_(0), noccB_(0), energy_(en), hcore_(make_shared<Hcore>(geom_,geom_->hcoreinfo())), nclosed_(_nclosed), nact_(_nact), nvirt_(_nvirt), nstate_(1), is_skelton_(false), ciwfn_(ci), rdm1_(_rdm1), rdm2_(_rdm2),
   rdm1_av_(_rdm1_av), rdm2_av_(_rdm2_av) {

  // we need to make sure that all the quantities are consistent in every MPI process
  if (c) {
    mpi__->broadcast(const_pointer_cast<Coeff>(c)->data(), c->size(), 0);
    coeff_ = c;
  }

  for (auto& i : *rdm1_)
    mpi__->broadcast(i.second->data(), i.second->size(), 0);
  for (auto& i : *rdm2_)
    mpi__->broadcast(i.second->data(), i.second->size(), 0);
  if (rdm1_av_)
    mpi__->broadcast(const_cast<double*>(rdm1_av_->data()), rdm1_av_->size(), 0);
  if (rdm2_av_)
    mpi__->broadcast(const_cast<double*>(rdm2_av_->data()), rdm2_av_->size(), 0);

  occup_ = VectorB(nclosed_+nact_+nvirt_);
  fill_n(occup_.data(), nclosed_, 2.0);

  //if (nact_ && rdm1_.empty())
  //  throw logic_error("If nact != 0, Reference::Reference wants to have RDMs.");

}


tuple<shared_ptr<const RDM<1>>,shared_ptr<const RDM<2>>> Reference::rdm12(const int ist, const int jst, const bool recompute) const {
  shared_ptr<const RDM<1>> r1;
  shared_ptr<const RDM<2>> r2;
  if (!recompute && rdm1_->exist(ist, jst) && rdm2_->exist(ist, jst)) {
    r1 = rdm1_->at(ist, jst);
    r2 = rdm2_->at(ist, jst);
  } else {
    FCI_bare fci(ciwfn_);
    fci.compute_rdm12(ist, jst);
    r1 = fci.rdm1(ist, jst);
    r2 = fci.rdm2(ist, jst);
  }
  return make_tuple(r1, r2);
}


tuple<shared_ptr<const RDM<3>>,shared_ptr<const RDM<4>>> Reference::rdm34(const int ist, const int jst) const {
  FCI_bare fci(ciwfn_);
  fci.compute_rdm12(ist, jst); // TODO stupid code
  return fci.rdm34(ist, jst);
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

shared_ptr<Matrix> Reference::rdm1_mat_tr(shared_ptr<const RDM<1>> active) const {
  if (nact_)
    return active->rdm1_mat_tr(nclosed_);
  else {
    auto out = make_shared<Matrix>(nocc(), nocc());
    for (int i = 0; i != nclosed_; ++i) out->element(i,i) = 0.0;	// unnecessary and silly: just to make it sure
    return out;
  }
}

shared_ptr<const Dvec> Reference::civectors() const {
  return ciwfn_->civectors();
}


shared_ptr<Dvec> Reference::rdm1deriv(const int istate) const {
  FCI_bare fci(ciwfn_);
  return fci.rdm1deriv(istate);
}


shared_ptr<Dvec> Reference::rdm2deriv(const int istate) const {
  FCI_bare fci(ciwfn_);
  return fci.rdm2deriv(istate);
}


shared_ptr<Matrix> Reference::rdm2fderiv(const int istate, shared_ptr<const Matrix> fock, shared_ptr<const Matrix> dmat) const {
  FCI_bare fci(ciwfn_);
  return fci.rdm2fderiv(istate, fock, dmat);
}

shared_ptr<Matrix> Reference::rdm2deriv_offset(const int istate, const size_t offset, const size_t size, shared_ptr<const Matrix> dmat) const {
  FCI_bare fci(ciwfn_);
  return fci.rdm2deriv_offset(istate, offset, size, dmat);
}


tuple<shared_ptr<Matrix>,shared_ptr<Matrix>>
Reference::rdm3deriv(const int istate, shared_ptr<const Matrix> fock, const size_t offset, const size_t size, shared_ptr<const Matrix> dbra_in, shared_ptr<const Matrix> fock_ebra_in) const {
  FCI_bare fci(ciwfn_);
  return fci.rdm3deriv(istate, fock, offset, size, dbra_in, fock_ebra_in);
}


shared_ptr<Reference> Reference::project_coeff(shared_ptr<const Geometry> geomin, const bool check_geom_change) const {

  if (geomin->magnetism())
    throw runtime_error("Projection from real to GIAO basis set is not implemented.   Use the GIAO code at zero-field.");

  bool moved = false;
  bool newbasis = false;

  if (check_geom_change) {
    auto j = geomin->atoms().begin();
    for (auto& i : geom_->atoms()) {
      moved |= i->distance(*j) > 1.0e-12;
      newbasis |= i->basis() != (*j)->basis();
      ++j;
    }
  } else {
    newbasis = true;
  }

  if (moved && newbasis)
    throw runtime_error("changing geometry and basis set at the same time is not allowed");

  shared_ptr<Reference> out;

  if (newbasis) {
    // project to a new basis
    const Overlap snew(geomin);
    Overlap snewinv = snew;
    snewinv.inverse_symmetric();
    MixedBasis<OverlapBatch> mixed(geom_, geomin);
    auto coeff = coeff_->copy();
    coeff->delocalize();
    auto c = make_shared<Coeff>(snewinv * mixed * *coeff);

    // make coefficient orthogonal (under the overlap metric)
    Matrix unit = *c % snew * *c;
    unit.inverse_half();
    *c *= unit;

    out = make_shared<Reference>(geomin, c, nclosed_, nact_, coeff_->mdim()-nclosed_-nact_, energy_);
    if (coeffA_) {
      assert(coeffB_);
      auto coeffA = coeffA_->copy();
      auto coeffB = coeffB_->copy();
      coeffA->delocalize();
      coeffB->delocalize();
      out->coeffA_ = make_shared<Coeff>(snewinv * mixed * *coeffA * unit);
      out->coeffB_ = make_shared<Coeff>(snewinv * mixed * *coeffB * unit);
    }
  } else {
    Overlap snew(geomin);
    Overlap sold(geom_);
    snew.inverse_half();
    sold.sqrt();
    auto coeff = coeff_->copy();
    coeff->delocalize();
    auto c = make_shared<Coeff>(snew * sold * *coeff);

    out = make_shared<Reference>(geomin, c, nclosed_, nact_, coeff_->mdim()-nclosed_-nact_, energy_);
    if (coeffA_) {
      assert(coeffB_);
      auto coeffA = coeffA_->copy();
      auto coeffB = coeffB_->copy();
      coeffA->delocalize();
      coeffB->delocalize();
      out->coeffA_ = make_shared<Coeff>(snew * sold * *coeffA);
      out->coeffB_ = make_shared<Coeff>(snew * sold * *coeffB);
    }
  }

  return out;
}


void Reference::set_eig(const VectorB& eig) {
  eig_ = eig;
  mpi__->broadcast(eig_.data(), eig_.size(), 0);
}


void Reference::set_eigB(const VectorB& eigB) {
  eigB_ = eigB;
  mpi__->broadcast(eigB_.data(), eigB_.size(), 0);
}


void Reference::set_occup(const VectorB& occup) {
  occup_ = occup;
  mpi__->broadcast(occup_.data(), occup_.size(), 0);
}


void Reference::set_occupB(const VectorB& occupB) {
  occupB_ = occupB;
  mpi__->broadcast(occupB_.data(), occupB_.size(), 0);
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


shared_ptr<Reference> Reference::extract_state(const vector<int> input, const bool update_rdms) const {
  cout << " * Extracting CI coefficients from Reference object for the following states: ";
  for (int i = 0; i != input.size(); ++i)
    cout << input[i] << " ";
  cout << endl;

  vector<double> newenergies(input.size());
  for (int i = 0; i != input.size(); ++i)
    newenergies[i] = energy_[input[i]];

  // Construct a CIWfn with only CI coefficients for the desired state
  auto newciwfn = make_shared<CIWfn>(geom_, nclosed_, nact_, input.size(), newenergies,
                                     ciwfn_->civectors()->extract_state(input), ciwfn_->det());

  // Use extract_average_rdm(...) to get desired RDMs and prepare output
  shared_ptr<Reference> out;
  if (update_rdms) {
    shared_ptr<Reference> rdmref = extract_average_rdm(input);
    out = make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, newenergies,
                                 rdmref->rdm1(), rdmref->rdm2(), rdmref->rdm1_av(), rdmref->rdm2_av(), newciwfn);
  } else {
    out = make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, newenergies,
                                 rdm1(), rdm2(), rdm1_av(), rdm2_av(), newciwfn);
  }
  return out;
}


// TODO Cleanup or remove?  Body is mostly the same as FCI_base::compute_rdm12()
shared_ptr<Reference> Reference::extract_average_rdm(const vector<int> rdm_state) const {
  if (rdm_state.size() == 0 || rdm_state.size() > nstate())
    throw runtime_error("Trying to obtain a state-averaged RDM over some invalid number of states.");

  cout << " * Extracting RDMs for ";
  cout << (rdm_state.size() > 1 ? "the average of the following states: " : "the following state: ");
  for (int i = 0; i != rdm_state.size(); ++i)
    cout << rdm_state[i] << " ";
  cout << endl;

  // for one-body RDM
  auto rdm1 = make_shared<VecRDM<1>>();
  auto rdm2 = make_shared<VecRDM<2>>();
  auto rdm1_av = make_shared<RDM<1>>(nact_);
  auto rdm2_av = make_shared<RDM<2>>(nact_);

  for (int index = 0; index != rdm_state.size(); ++index) {
    const int istate = rdm_state[index];
    // one and two body RDMs
    rdm1->emplace(index, rdm1_->at(istate)->copy());
    rdm2->emplace(index, rdm2_->at(istate)->copy());
  }

  if (rdm_state.size() > 1) {
    for (int index = 0; index != rdm_state.size(); ++index) {
      const int istate = rdm_state[index];
      const double weight = 1.0/static_cast<double>(rdm_state.size());
      rdm1_av->ax_plus_y(weight, rdm1_->at(istate));
      rdm2_av->ax_plus_y(weight, rdm2_->at(istate));
    }
  } else {
    rdm1_av = rdm1_->at(rdm_state[0])->copy();
    rdm2_av = rdm2_->at(rdm_state[0])->copy();
  }

  return make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, energy_, rdm1, rdm2, rdm1_av, rdm2_av, ciwfn_);
}


