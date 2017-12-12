//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/rasci.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/ci/fci/modelci.h>
#include <src/ci/ras/rasci.h>
#include <src/ci/ras/form_sigma.h>
#include <src/util/combination.hpp>
#include <src/util/math/davidson.h>

using namespace std;
using namespace bagel;

RASCI::RASCI(shared_ptr<const PTree> idat, shared_ptr<const Geometry> g, shared_ptr<const Reference> r)
 : Method(idat, g, r) {
  common_init();
  update(ref_->coeff());
}

void RASCI::common_init() {
  print_header();

//const bool frozen = idata_->get<bool>("frozen", false);
  max_iter_ = idata_->get<int>("maxiter", 100);
  davidson_subspace_ = idata_->get<int>("davidson_subspace", 20);
  thresh_ = idata_->get<double>("thresh", 1.0e-8);
  print_thresh_ = idata_->get<double>("print_thresh", 0.05);

  batchsize_ = idata_->get<int>("batchsize", 512);

  nstate_ = idata_->get<int>("nstate", 1);
  nguess_ = idata_->get<int>("nguess", nstate_);

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
  const bool extern_nactele = idata_->get<bool>("extern_nactele", false);
  if (!extern_nactele) {
    if ((geom_->nele()+nspin-charge) % 2 != 0) throw runtime_error("Invalid nspin specified in RASCI");
    nelea_ = (geom_->nele()+nspin-charge)/2 - ncore_;
    neleb_ = (geom_->nele()-nspin-charge)/2 - ncore_;
  } else {
    const int nactele = idata_->get<int>("nactele");
    nelea_ = (nactele+nspin-charge) / 2;
    if ((nactele+nspin-charge) % 2 != 0) throw runtime_error("Invalid nspin specified in RASCI");
    neleb_ = nactele - charge - nelea_;
    assert(neleb_ == (nactele-nspin-charge)/2);
  }

  // TODO allow for zero electron (quick return)
  if (nelea_ < 0 || neleb_ < 0) throw runtime_error("#electrons cannot be negative in RASCI");
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

void RASCI::model_guess(shared_ptr<RASDvec>& out) {
  multimap<double, pair<bitset<nbit__>, bitset<nbit__>>> ordered_elements;
  for (auto& b : denom_->blocks()) {
    if (!b) continue;
    const double* d = b->data();
    for (auto& abit : *b->stringsa()) {
      for (auto& bbit : *b->stringsb()) {
        ordered_elements.emplace(*d++, make_pair(abit, bbit));
      }
    }
  }

  vector<pair<bitset<nbit__>, bitset<nbit__>>> basis;
  double last_value = 0.0;
  for (auto& p : ordered_elements) {
    double val = p.first;
    if (basis.size() >= nguess_ && val != last_value)
      break;
    else
      basis.push_back(p.second);
  }
  const int nguess = basis.size();

  shared_ptr<Matrix> spin = make_shared<CISpin>(basis, norb_);
  VectorB eigs(nguess);
  spin->diagonalize(eigs);

  int start, end;
  const double target_spin = 0.25 * static_cast<double>(det_->nspin()*(det_->nspin()+2));
  for (start = 0; start < nguess; ++start)
    if (fabs(eigs(start) - target_spin) < 1.0e-8) break;
  for (end = start; end < nguess; ++end)
    if (fabs(eigs(end) - target_spin) > 1.0e-8) break;

  if ((end-start) >= nstate_) {
    const MatView coeffs = spin->slice(start, end);

    shared_ptr<Matrix> hamiltonian = make_shared<CIHamiltonian>(basis, jop_);
    hamiltonian = make_shared<Matrix>(coeffs % *hamiltonian * coeffs);
    hamiltonian->diagonalize(eigs);

#if 0
    const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();
    for (int i = 0; i < end-start; ++i)
      cout << setw(12) << setprecision(8) << eigs(i) + nuc_core << endl;
#endif

    auto coeffs1 = make_shared<Matrix>(coeffs * *hamiltonian);
    for (int i = 0; i < nguess; ++i) {
      const bitset<nbit__> ia = basis[i].first;
      const bitset<nbit__> ib = basis[i].second;

      for (int j = 0; j < nstate_; ++j)
        out->data(j)->element(ib, ia) = coeffs1->element(i, j);
    }
  }
  else if (nguess_ >= det_->size()) {
    stringstream message;
    message << "Asking for " << nstate_ << " states, but there seems to only be " << end-start << " states with the right spin.";
    throw runtime_error(message.str());
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

    // This can happen if all possible determinants are checked without finding nstate acceptable ones.
    if (alpha.count() + beta.count() != nelea_ + neleb_)
      throw logic_error("RasCI::generate_guess produced an invalid determinant.  Check the number of states being requested.");

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
          setw(20) << left << print_bit(alpha&beta, det()->norb()) << " open " << setw(20) << print_bit(open_bit, det()->norb()) << right << endl;

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
  for (int i = 0; i != ndet; ++i) tmp.emplace(-1.0e10*(1+i), make_pair(bitset<nbit__>(0),bitset<nbit__>(0)));

  for (auto& iblock : denom_->blocks()) {
    if (!iblock) continue;
    double* diter = iblock->data();
    for (auto& aiter : *iblock->stringsa()) {
      for (auto& biter : *iblock->stringsb()) {
        const double din = -(*diter);
        if (tmp.begin()->first < din) {
          tmp.emplace(din, make_pair(biter, aiter));
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

  // Creating an initial CI vector
  cc_ = make_shared<RASDvec>(det_, nstate_);

  // find determinants that have small diagonal energies
  if (nguess_ <= nstate_)
    generate_guess(nelea_-neleb_, nstate_, cc_);
  else
    model_guess(cc_);
  pdebug.tick_print("guess generation");

  // nuclear energy retrieved from geometry
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();

  // Davidson utility
  DavidsonDiag<RASCivec> davidson(nstate_, davidson_subspace_);

  // Object in charge of forming sigma vector
  FormSigmaRAS form_sigma(batchsize_);

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
      else {
        ccn.push_back(nullptr);
        sigman.push_back(nullptr);
      }
    }
    const vector<double> energies = davidson.compute(ccn, sigman);

    // get residual and new vectors
    vector<shared_ptr<RASCivec>> errvec = davidson.residual();
    pdebug.tick_print("davidson");

    // compute errors
    vector<double> errors;
    for (int i = 0; i != nstate_; ++i) {
      errors.push_back(errvec[i]->rms());
      conv[i] = static_cast<int>(errors[i] < thresh_);
    }
    pdebug.tick_print("error");

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
        cc_->data(ist)->normalize();
        cc_->data(ist)->spin_decontaminate();
        cc_->data(ist)->synchronize();
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
