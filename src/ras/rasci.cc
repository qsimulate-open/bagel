//
// BAGEL - Parallel electron correlation program.
// Filename: ras/rasci.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/ras/rasci.h>
#include <src/ras/form_sigma.h>
#include <src/util/combination.hpp>
#include <src/math/davidson.h>

using namespace std;
using namespace bagel;

RASCI::RASCI(shared_ptr<const PTree> idat, shared_ptr<const Geometry> g, shared_ptr<const Reference> r)
 : Method(idat, g, r) {
  common_init();
  update(ref_->coeff());
}

void RASCI::common_init() {
  print_header();

  const bool frozen = idata_->get<bool>("frozen", false);
  max_iter_ = idata_->get<int>("maxiter", 100);
  davidsonceiling_ = idata_->get<int>("davidsonceiling", 10);
  thresh_ = idata_->get<double>("thresh", 1.0e-16);
  print_thresh_ = idata_->get<double>("print_thresh", 0.05);
  sparse_ = idata_->get<bool>("sparse", true);

  nstate_ = idata_->get<int>("nstate", 1);

  // No defaults for RAS, must set "active"
  const shared_ptr<const PTree> iactive = idata_->get_child("active");
  if (iactive->size() != 3) throw runtime_error("Must specify three active spaces in RAS calculations.");
  vector<set<int>> acts;
  for (auto& i : *iactive) {
    set<int> tmpset;
    for (auto& j : *i)
      if (!tmpset.insert(lexical_cast<int>(j->data()) - 1).second) throw runtime_error("Duplicate orbital in list of active orbitals.");
    acts.push_back(tmpset);
  }
  ref_ = ref_->set_ractive(acts[0], acts[1], acts[2]);
  ncore_ = ref_->nclosed();

  ras_ = {{ static_cast<int>(acts[0].size()), static_cast<int>(acts[1].size()), static_cast<int>(acts[2].size()) }};
  norb_ = ras_[0] + ras_[1] + ras_[2];

  max_holes_ = idata_->get<int>("max_holes", 0);
  max_particles_ = idata_->get<int>("max_particles", 0);

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
  if (nelea_ <= 0 || neleb_ <= 0) throw runtime_error("#electrons cannot be zero/negative in RASCI");
  //for (int i = 0; i != nstate_; ++i) weight_.push_back(1.0/static_cast<double>(nstate_));

#ifndef NORDMS
  // resizing rdm vectors (with null pointers)
  rdm1_.resize(nstate_);
  rdm2_.resize(nstate_);
#endif
  energy_.resize(nstate_);

  // construct a determinant space in which this RASCI will be performed.
  det_ = make_shared<const RASDeterminants>(ras_, nelea_, neleb_, max_holes_, max_particles_);
}

// generate initial vectors
//   - bits: bit patterns of low-energy determinants
//   - nspin: #alpha - #beta
//   - out:
void RASCI::generate_guess(const int nspin, const int nstate, shared_ptr<RASDvec>& out) {
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

    pair<vector<tuple<bitset<nbit__>, bitset<nbit__>, int>>, double> adapt = det()->spin_adapt(nelea_-neleb_, alpha, beta);
    const double fac = adapt.second;
    for (auto& iter : adapt.first) {
      out->data(oindex)->element(get<0>(iter), get<1>(iter)) = get<2>(iter)*fac;
    }
    out->data(oindex)->spin_decontaminate();

    cout << "     guess " << setw(3) << oindex << ":   closed " <<
          setw(20) << left << det()->print_bit(alpha&beta) << " open " << setw(20) << det()->print_bit(open_bit) << right << endl;

    ++oindex;
    if (oindex == nstate) break;
  }
  if (oindex < nstate) {
    for (auto& io : out->dvec()) io->zero();
    ndet *= 4;
    goto start_over;
  }
  cout << endl;
}

// returns seed determinants for initial guess
vector<pair<bitset<nbit__> , bitset<nbit__>>> RASCI::detseeds(const int ndet) {
  multimap<double, pair<bitset<nbit__>,bitset<nbit__>>> tmp;
  for (int i = 0; i != ndet; ++i) tmp.insert(make_pair(-1.0e10*(1+i), make_pair(bitset<nbit__>(0),bitset<nbit__>(0))));

  for (auto& iblock : denom_->blocks()) {
    if (!iblock) continue;
    double* diter = iblock->data();
    for (auto& aiter : *iblock->stringa()) {
      for (auto& biter : *iblock->stringb()) {
        const double din = -(*diter);
        if (tmp.begin()->first < din) {
          tmp.insert(make_pair(din, make_pair(biter, aiter)));
          tmp.erase(tmp.begin());
        }
        ++diter;
      }
    }
  }
  assert(tmp.size() == ndet || ndet > det_->size());
  vector<pair<bitset<nbit__> , bitset<nbit__>>> out;
  for (auto iter = tmp.rbegin(); iter != tmp.rend(); ++iter)
    out.push_back(iter->second);
  return out;
}

void RASCI::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "        RASCI calculation      " << endl;
  cout << "  ---------------------------" << endl << endl;
}


void RASCI::compute() {
  Timer pdebug(0);

  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("RASCI: C1 only at the moment.");

  // Creating an initial CI vector
  cc_ = make_shared<RASDvec>(det_, nstate_);

  // find determinants that have small diagonal energies
  generate_guess(nelea_-neleb_, nstate_, cc_);
  pdebug.tick_print("guess generation");

  // nuclear energy retrieved from geometry
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();

  // Davidson utility
  DavidsonDiag<RASCivec> davidson(nstate_, davidsonceiling_);

  // Object in charge of forming sigma vector
  FormSigmaRAS form_sigma(sparse_);

  // main iteration starts here
  cout << "  === RAS-CI iteration ===" << endl << endl;
  // 0 means not converged
  vector<int> conv(nstate_,0);

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer fcitime;

    // form a sigma vector given cc
    shared_ptr<RASDvec> sigma = form_sigma(cc_, jop_, conv);
    pdebug.tick_print("sigma vector");

    // constructing Dvec's for Davidson
    vector<shared_ptr<const RASCivec>> ccn, sigman;
    for (int i = 0; i < nstate_; ++i) {
      if (!conv[i]) {
        ccn.push_back(make_shared<const RASCivec>(*cc_->data(i)));
        sigman.push_back(make_shared<const RASCivec>(*sigma->data(i)));
      }
    }
    const vector<double> energies = davidson.compute(ccn, sigman);

    // get residual and new vectors
    vector<shared_ptr<RASCivec>> errvec = davidson.residual();
    pdebug.tick_print("davidson");

    // compute errors
    vector<double> errors;
    for (int i = 0; i != nstate_; ++i) {
      errors.push_back(errvec[i]->variance());
      conv[i] = static_cast<int>(errors[i] < thresh_);
    }
    pdebug.tick_print("error");

#ifdef HAVE_MPI_H
    vector<int> requests;
#endif

    if (!*min_element(conv.begin(), conv.end())) {
      // denominator scaling
      for (int ist = 0; ist != nstate_; ++ist) {
        if (conv[ist]) continue;
        const size_t size = cc_->data(ist)->size();
        double* target_array = cc_->data(ist)->data();
        double* source_array = errvec.at(ist)->data();
        double* denom_array = denom_->data();
        const double en = energies.at(ist);
        transform(source_array, source_array + size, denom_array, target_array, [&en] (const double cc, const double den) { return cc / min(en - den, -0.1); });
        davidson.orthog(cc_->data(ist));
        list<shared_ptr<const RASCivec>> tmp;
        for (int jst = 0; jst != ist; ++jst) tmp.push_back(cc_->data(jst));
        cc_->data(ist)->orthog(tmp);
        cc_->data(ist)->spin_decontaminate();
#ifdef HAVE_MPI_H
        requests.push_back(mpi__->ibroadcast(cc_->data(ist)->data(), cc_->data(ist)->size(), 0));
#endif
      }
    }
#ifdef HAVE_MPI_H
    for (auto& r : requests)
      mpi__->wait(r);
#endif

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

  cc_ = make_shared<RASDvec>(davidson.civec());

  for (int istate = 0; istate < nstate_; ++istate) {
    cout << endl << "     * ci vector " << setw(3) << istate << ", <S^2> = " << setw(6) << setprecision(4) << cc_->data(istate)->spin_expectation()
                 << ", E = " << setw(17) << fixed << setprecision(8) << energy_[istate] << endl;
    cc_->data(istate)->print(print_thresh_);
  }

#if 0
  for (auto& iprop : properties_) {
    iprop->compute(cc_);
    iprop->print();
  }
#endif
}
