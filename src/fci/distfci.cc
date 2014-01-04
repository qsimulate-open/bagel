//
// BAGEL - Parallel electron correlation program.
// Filename: fci/distfci.cc
// Copyright (C) 2011 Toru Shiozaki
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

#include <cassert>

#include <src/util/combination.hpp>
#include <src/math/comb.h>
#include <src/fci/modelci.h>
#include <src/fci/dist_form_sigma.h>
#include <src/fci/distfci.h>
#include <src/math/davidson.h>
#include <src/fci/space.h>
#include <src/fci/hzdenomtask.h>

using namespace std;
using namespace bagel;

DistFCI::DistFCI(std::shared_ptr<const PTree> idat, shared_ptr<const Geometry> g, shared_ptr<const Reference> r, const int ncore, const int norb, const int nstate)
 : Method(idat, g, r), ncore_(ncore), norb_(norb), nstate_(nstate) {
  common_init();

#ifndef HAVE_MPI_H
  throw logic_error("DistFCI can be used only with MPI");
#endif

  cout << "    * Parallel algorithm will be used." << endl;

  update(ref_->coeff());
}

void DistFCI::common_init() {
  print_header();

  const bool frozen = idata_->get<bool>("frozen", false);
  max_iter_ = idata_->get<int>("maxiter", 100);
  max_iter_ = idata_->get<int>("maxiter_fci", max_iter_);
  davidsonceiling_ = idata_->get<int>("davidsonceiling", 10);
  thresh_ = idata_->get<double>("thresh", 1.0e-20);
  thresh_ = idata_->get<double>("thresh_fci", thresh_);
  print_thresh_ = idata_->get<double>("print_thresh", 0.05);

  if (nstate_ < 0) nstate_ = idata_->get<int>("nstate", 1);
  nguess_ = idata_->get<int>("nguess", nstate_);

  const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
  if (iactive) {
    set<int> tmp;
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iactive) tmp.insert(lexical_cast<int>(i->data()) - 1);
    ref_ = ref_->set_active(tmp);
    ncore_ = ref_->nclosed();
    norb_ = ref_->nact();
  }
  else {
    if (ncore_ < 0) ncore_ = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
    if (norb_  < 0) norb_ = idata_->get<int>("norb", ref_->coeff()->mdim()-ncore_);
  }

  // Configure properties to be calculated on the final wavefunctions
  //if (idata_->get<bool>("dipoles", false)) properties_.push_back(make_shared<CIDipole>(ref_, ncore_, ncore_+norb_));

  // additional charge
  const int charge = idata_->get<int>("charge", 0);

  // nspin is #unpaired electron 0:singlet, 1:doublet, 2:triplet, ... (i.e., Molpro convention).
  const int nspin = idata_->get<int>("nspin", 0);
  if ((geom_->nele()+nspin-charge) % 2 != 0) throw runtime_error("Invalid nspin specified");
  nelea_ = (geom_->nele()+nspin-charge)/2 - ncore_;
  neleb_ = (geom_->nele()-nspin-charge)/2 - ncore_;

  // TODO allow for zero electron (quick return)
  if (nelea_ <= 0 || neleb_ <= 0) throw runtime_error("#electrons cannot be zero/negative in FCI");
  //for (int i = 0; i != nstate_; ++i) weight_.push_back(1.0/static_cast<double>(nstate_));

  // resizing rdm vectors (with null pointers)
  //rdm1_.resize(nstate_);
  //rdm2_.resize(nstate_);
  energy_.resize(nstate_);

  // construct a determinant space in which this FCI will be performed.
  space_ = make_shared<Space>(norb_, nelea_, neleb_, 1);
  det_ = space_->basedet();
}

void DistFCI::model_guess(vector<shared_ptr<DistCivec>>& out) {
  multimap<double, pair<size_t, size_t>> ordered_elements;
  const double* d = denom_->local();
  for (size_t ia = denom_->astart(); ia < denom_->aend(); ++ia) {
    for (size_t ib = 0; ib < det_->lenb(); ++ib) {
      ordered_elements.emplace(*d++, make_pair(ia, ib));
    }
  }

  vector<double> energies;
  vector<size_t> aarray, barray;
  double last_value = 0.0;
  for (auto& p : ordered_elements) {
    double val = p.first;
    if (energies.size() >= nguess_ && val != last_value)
      break;
    else {
      energies.push_back(val);
      aarray.push_back(p.second.first);
      barray.push_back(p.second.second);
    }
  }

  vector<size_t> nelements(mpi__->size(), 0);
  const size_t nn = energies.size();
  mpi__->allgather(&nn, 1, nelements.data(), 1);

  const size_t chunk = *max_element(nelements.begin(), nelements.end());
  energies.resize(chunk, 0);
  aarray.resize(chunk, 0);
  barray.resize(chunk, 0);

  vector<double> allenergies(chunk * mpi__->size(), 0.0);
  mpi__->allgather(energies.data(), chunk, allenergies.data(), chunk);
  vector<size_t> allalpha(chunk * mpi__->size());
  mpi__->allgather(aarray.data(), chunk, allalpha.data(), chunk);
  vector<size_t> allbeta(chunk * mpi__->size());
  mpi__->allgather(barray.data(), chunk, allbeta.data(), chunk);

  ordered_elements.clear();
  for (size_t i = 0; i < chunk * mpi__->size(); ++i) {
    if (allenergies[i] != 0.0) ordered_elements.emplace(allenergies[i], make_pair(allalpha[i], allbeta[i]));
  }

  vector<pair<bitset<nbit__>, bitset<nbit__>>> basis;
  last_value = 0.0;
  for (auto& p : ordered_elements) {
    double val = p.first;
    if (basis.size() >= nguess_ && val != last_value)
      break;
    else
      basis.emplace_back(det_->stringa(p.second.first), det_->stringb(p.second.second));
  }
  const int nguess = basis.size();

  shared_ptr<Matrix> spin = make_shared<CISpin>(basis, norb_);
  vector<double> eigs(nguess, 0.0);
  spin->diagonalize(eigs.data());

  int start, end;
  const double target_spin = 0.25 * static_cast<double>(det_->nspin()*(det_->nspin()+2));
  for (start = 0; start < nguess; ++start)
    if (fabs(eigs[start] - target_spin) < 1.0e-8) break;
  for (end = start; end < nguess; ++end)
    if (fabs(eigs[end] - target_spin) > 1.0e-8) break;

  if ((end-start) >= nstate_) {
    shared_ptr<Matrix> coeffs = spin->slice(start, end);

    shared_ptr<Matrix> hamiltonian = make_shared<CIHamiltonian>(basis, jop_);
    hamiltonian = make_shared<Matrix>(*coeffs % *hamiltonian * *coeffs);
    hamiltonian->diagonalize(eigs.data());

#if 0
    const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();
    for (int i = 0; i < end-start; ++i)
      cout << setw(12) << setprecision(8) << eigs[i] + nuc_core << endl;
#endif

    coeffs = (*coeffs * *hamiltonian).slice(0, nstate_);
    mpi__->broadcast(coeffs->data(), coeffs->ndim() * coeffs->mdim(), 0);
    const size_t lenb = det_->lenb();
    for (int i = 0; i < nguess; ++i) {
      size_t ia = det_->lexical<0>(basis[i].first);
      if ( ia >= denom_->astart() && ia < denom_->aend() ) {
        ia -= denom_->astart();
        const size_t ib = det_->lexical<1>(basis[i].second);
        for (int j = 0; j < nstate_; ++j)
          out[j]->local(ib + ia*lenb) = coeffs->element(i, j);
      }
    }
  }
  else {
    nguess_ *= 2;
    model_guess(out);
  }
}

// generate initial vectors
//   - bits: bit patterns of low-energy determinants
//   - nspin: #alpha - #beta
//   - out:
void DistFCI::generate_guess(const int nspin, const int nstate, vector<shared_ptr<DistCivec>>& out) {
  int ndet = nstate_*10;
  start_over:
  vector<pair<bitset<nbit__>, bitset<nbit__>>> bits = detseeds(ndet);

  // Spin adapt detseeds
  int oindex = 0;
  vector<bitset<nbit__>> done;
  for (auto& it : bits) {
    bitset<nbit__> alpha = it.second;
    bitset<nbit__> beta = it.first;
    bitset<nbit__> open_bit = (alpha^beta);

    // make sure that we have enough unpaired alpha
    const int unpairalpha = (alpha ^ (alpha & beta)).count();
    const int unpairbeta  = (beta ^ (alpha & beta)).count();
    if (unpairalpha-unpairbeta < nelea_-neleb_) continue;

    // check if this orbital configuration is already used
    if (find(done.begin(), done.end(), open_bit) != done.end()) continue;
    done.push_back(open_bit);

    pair<vector<tuple<int, int, int>>, double> adapt = det()->spin_adapt(nelea_-neleb_, alpha, beta);
    const double fac = adapt.second;

    out[oindex]->zero();
    for (auto& ad : adapt.first) {
      const int aloc = get<1>(ad) - out[oindex]->astart();
      if (aloc >= 0 && aloc < out[oindex]->asize())
        out[oindex]->local(get<0>(ad) + det_->lenb()*aloc) = get<2>(ad)*fac;
    }
    out[oindex]->spin_decontaminate();

    cout << "     guess " << setw(3) << oindex << ":   closed " <<
          setw(20) << left << det()->print_bit(alpha&beta) << " open " << setw(20) << det()->print_bit(open_bit) << right << endl;

    ++oindex;
    if (oindex == nstate) break;
  }
  if (oindex < nstate) {
    for (auto& i : out) i->zero();
    ndet *= 4;
    goto start_over;
  }
  cout << endl;
}

vector<pair<bitset<nbit__> , bitset<nbit__>>> DistFCI::detseeds(const int ndet) {
  multimap<double, pair<size_t, size_t>> tmp;
  for (int i = 0; i != ndet; ++i)
    tmp.insert(make_pair(-1.0e10*(1+i), make_pair(0,0)));

  double* diter = denom_->local();
  for (size_t ia = denom_->astart(); ia != denom_->aend(); ++ia) {
    for (size_t ib = 0; ib != det_->lenb(); ++ib) {
      const double din = -*diter++;
      if (tmp.begin()->first < din) {
        tmp.insert(make_pair(din, make_pair(ib, ia)));
        tmp.erase(tmp.begin());
      }
    }
  }

  assert(ndet == tmp.size());

  vector<size_t> aarray, barray;
  vector<double> en;
  for (auto iter = tmp.rbegin(); iter != tmp.rend(); ++iter) {
    aarray.push_back(iter->second.second);
    barray.push_back(iter->second.first);
    en.push_back(iter->first);
  }

  // rank 0 will take care of this
  vector<size_t> aall(mpi__->size()*ndet);
  vector<size_t> ball(mpi__->size()*ndet);
  vector<double> eall(mpi__->size()*ndet);
  mpi__->allgather(aarray.data(), ndet, aall.data(), ndet);
  mpi__->allgather(barray.data(), ndet, ball.data(), ndet);
  mpi__->allgather(en.data(),     ndet, eall.data(), ndet);

  tmp.clear();
  for (int i = 0; i != aall.size(); ++i) {
    tmp.insert(make_pair(eall[i], make_pair(ball[i], aall[i])));
  }

  // sync'ing to make sure the consistency
  auto c = tmp.rbegin();
  for (int i = 0; i != ndet; ++i, ++c) {
    ball[i] = c->second.first;
    aall[i] = c->second.second;
  }
  mpi__->broadcast(aall.data(), ndet, 0);
  mpi__->broadcast(ball.data(), ndet, 0);

  vector<pair<bitset<nbit__> , bitset<nbit__>>> out;
  for (int i = 0; i != ndet; ++i)
    out.push_back(make_pair(det_->stringb(ball[i]), det_->stringa(aall[i])));

  return out;
}


void DistFCI::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "        FCI calculation      " << endl;
  cout << "  ---------------------------" << endl << endl;
}



void DistFCI::update(shared_ptr<const Coeff> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  jop_ = make_shared<Jop>(ref_, ncore_, ncore_+norb_, c, "HZ");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;

  const_denom();
}



// same as HZ::const_denom except that denom_ is also distributed
void DistFCI::const_denom() {
  Timer denom_t;
  auto h = make_shared<Matrix>(norb_, 1);
  auto jop = make_shared<Matrix>(norb_, norb_);
  auto kop = make_shared<Matrix>(norb_, norb_);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop->element(i,j) = jop->element(j,i) = 0.5*jop_->mo2e_hz(i, j, i, j);
      kop->element(i,j) = kop->element(j,i) = 0.5*jop_->mo2e_hz(i, j, j, i);
    }
    h->element(i,0) = jop_->mo1e(i,i);
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<DistCivec>(det_);

  double* iter = denom_->local();
  TaskQueue<HZDenomTask> tasks(denom_->asize());
  for (size_t i = denom_->astart(); i != denom_->aend(); ++i) {
    tasks.emplace_back(iter, denom_->det()->stringa(i), det_, jop, kop, h);
    iter += det()->stringb().size();
  }
  tasks.compute();

  denom_t.tick_print("denom");
}

void DistFCI::compute() {
  Timer pdebug(0);

  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("FCI: C1 only at the moment.");

  // some constants
  //const int ij = nij();

  // Creating an initial CI vector
  vector<shared_ptr<DistCivec>> cc(nstate_);
  for (auto& i : cc)
    i = make_shared<DistCivec>(det_);

  // find determinants that have small diagonal energies
  if (nguess_ <= nstate_)
    generate_guess(nelea_-neleb_, nstate_, cc);
  else
    model_guess(cc);
  pdebug.tick_print("guess generation");

  // nuclear energy retrieved from geometry
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();

  // Davidson utility
  DavidsonDiag<DistCivec> davidson(nstate_, davidsonceiling_);

  // main iteration starts here
  cout << "  === FCI iteration ===" << endl << endl;
  // 0 means not converged
  vector<int> conv(nstate_, 0);

  FormSigmaDistFCI form_sigma(space_);

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer fcitime;

    // form a sigma vector given cc
    vector<shared_ptr<DistCivec>> sigma = form_sigma(cc, jop_, conv);
    pdebug.tick_print("sigma vector");

    // Davidson
    vector<shared_ptr<const DistCivec>> ccn, sigman;
    for (auto& i : cc) if (i) ccn.push_back(i);
    for (auto& i : sigma) if (i) sigman.push_back(i);
    const vector<double> energies = davidson.compute(ccn, sigman);

    // get residual and new vectors
    vector<shared_ptr<DistCivec>> errvec = davidson.residual();
    pdebug.tick_print("davidson");

    // compute errors
    vector<double> errors;
    for (int i = 0; i != nstate_; ++i) {
      errors.push_back(errvec[i]->variance());
      conv[i] = static_cast<int>(errors[i] < thresh_);
    }
    pdebug.tick_print("error");

    cc.clear();
    if (!*min_element(conv.begin(), conv.end())) {
      // denominator scaling
      for (int ist = 0; ist != nstate_; ++ist) {
        if (!conv[ist]) {
          shared_ptr<DistCivec> c = errvec[ist]->clone();
          const int size = c->size();
          double* target_array = c->local();
          double* source_array = errvec[ist]->local();
          double* denom_array = denom_->local();
          const double en = energies[ist];
          // TODO this should be threaded
          for (int i = 0; i != size; ++i) {
            target_array[i] = source_array[i] / min(en - denom_array[i], -0.1);
          }
          c->spin_decontaminate();
          davidson.orthog(c);
          list<shared_ptr<const DistCivec>> tmp;
          for (int jst = 0; jst != ist; ++jst)
            if (!conv[jst]) tmp.push_back(cc.at(jst));
          c->orthog(tmp);
          cc.push_back(c);
        } else {
          cc.push_back(shared_ptr<DistCivec>());
        }
      }
    }
    pdebug.tick_print("denominator");

    // printing out
    if (nstate_ != 1 && iter) cout << endl;
    for (int i = 0; i != nstate_; ++i) {
      cout << setw(7) << iter << setw(3) << i << setw(2) << (conv[i] ? "*" : " ")
                              << setw(17) << fixed << setprecision(8) << energies[i]+nuc_core << "   "
                              << setw(10) << scientific << setprecision(2) << errors[i] << fixed << setw(10) << setprecision(2)
                              << fcitime.tick() << endl;
      energy_[i] = energies[i]+nuc_core;
    }
    if (*min_element(conv.begin(), conv.end())) break;
  }
  // main iteration ends here

  // TODO RDM etc is not properly done yet
  cc_ = make_shared<DistDvec>(davidson.civec());
  for (int ist = 0; ist < nstate_; ++ist) {
    const double s2 = cc_->data(ist)->spin_expectation();
    if (mpi__->rank() == 0)
      cout << endl << "     * ci vector " << setw(3) << ist
                   << ", <S^2> = " << setw(6) << setprecision(4) << s2
                   << ", E = " << setw(17) << fixed << setprecision(8) << energy_[ist] << endl;
    cc_->data(ist)->print(print_thresh_);
  }
}
