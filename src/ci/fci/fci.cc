//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci.cc
// Copyright (C) 2011 Toru Shiozaki
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

#include <src/ci/fci/fci.h>
#include <src/ci/fci/harrison.h>
#include <src/ci/fci/knowles.h>
#include <src/ci/fci/space.h>
#include <src/ci/fci/modelci.h>
#include <src/util/combination.hpp>
#include <src/util/exception.h>
#include <src/prop/multipole.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(FCI)

FCI::FCI(shared_ptr<const PTree> idat, shared_ptr<const Geometry> g, shared_ptr<const Reference> r,
         const int ncore, const int norb, const int nstate, const bool store)
 : FCI_base(idat, g, r, ncore, norb, nstate, store) {
  common_init();
}

HarrisonZarrabian::HarrisonZarrabian(shared_ptr<const PTree> a, shared_ptr<const Geometry> g, shared_ptr<const Reference> b,
                                     const int ncore, const int nocc, const int nstate, const bool store) : FCI(a, g, b, ncore, nocc, nstate, store) {
  space_ = make_shared<HZSpace>(det_);
  update(ref_->coeff());
  if (idata_->get<bool>("only_ints", false)) {
    OArchive ar("ref");
    ar << ref_;
    dump_ints();
    throw Termination("MO integrals are dumped on a file.");
  }
}


KnowlesHandy::KnowlesHandy(shared_ptr<const PTree> a, shared_ptr<const Geometry> g, shared_ptr<const Reference> b,
                           const int ncore, const int nocc, const int nstate, const bool store) : FCI(a, g, b, ncore, nocc, nstate, store) {
  update(ref_->coeff());
  if (idata_->get<bool>("only_ints", false)) {
    OArchive ar("ref");
    ar << ref_;
    dump_ints();
    throw Termination("MO integrals are dumped on a file.");
  }
}


KnowlesHandy::KnowlesHandy(shared_ptr<const CIWfn> ci, shared_ptr<const Reference> r) {
  print_thresh_ = 1.0e-8;
  nelea_ = ci->det()->nelea();
  neleb_ = ci->det()->neleb();
  ncore_ = ci->ncore();
  norb_  = ci->nact();
  nstate_ = ci->nstates();
  energy_ = ci->energies();
  cc_ = ci->civectors()->copy();
// Since the determinant space might not be compatible, reconstruct determinant space
  det_ = make_shared<const Determinants>(norb_, nelea_, neleb_, /*compress=*/true, /*mute=*/true);
  rdm1_ = make_shared<VecRDM<1>>();
  rdm2_ = make_shared<VecRDM<2>>();
  ref_ = r;
  store_half_ints_ = false;
  weight_ = vector<double>(nstate_, 1.0/static_cast<double>(nstate_));
  update(ref_->coeff());
}



void FCI::common_init() {
  print_header();

  const bool frozen = idata_->get<bool>("frozen", false);
  max_iter_ = idata_->get<int>("maxiter", 100);
  max_iter_ = idata_->get<int>("maxiter_fci", max_iter_);
  davidson_subspace_ = idata_->get<int>("davidson_subspace", 20);
  thresh_ = idata_->get<double>("thresh", 1.0e-10);
  thresh_ = idata_->get<double>("thresh_fci", thresh_);
  print_thresh_ = idata_->get<double>("print_thresh", 0.05);
  restart_ = idata_->get<bool>("restart", false);

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

  // calculate dipole moments if requested
  dipoles_ = idata_->get<bool>("dipoles", false);

  // additional charge
  const int charge = idata_->get<int>("charge", 0);

  // nspin is #unpaired electron 0:singlet, 1:doublet, 2:triplet, ... (i.e., Molpro convention).
  const int nspin = idata_->get<int>("nspin", 0);
  if ((geom_->nele()+nspin-charge) % 2 != 0)
    throw runtime_error("Invalid nspin specified");
  nelea_ = (geom_->nele()+nspin-charge)/2 - ncore_;
  neleb_ = (geom_->nele()-nspin-charge)/2 - ncore_;

  // TODO allow for zero electron (quick return)
  // neleb can be = 0, so long as nelea > 0
  if (nelea_ <= 0 || neleb_ < 0)
    throw runtime_error("#electrons cannot be zero/negative in FCI");
  weight_ = vector<double>(nstate_, 1.0/static_cast<double>(nstate_));

  // initialize VecRDM
  rdm1_ = make_shared<VecRDM<1>>();
  rdm2_ = make_shared<VecRDM<2>>();
  energy_.resize(nstate_);

  // construct a determinant space in which this FCI will be performed.
  det_ = make_shared<const Determinants>(norb_, nelea_, neleb_);

}


void FCI::model_guess(shared_ptr<Dvec> out) {
  multimap<double, pair<bitset<nbit__>, bitset<nbit__>>> ordered_elements;
  const double* d = denom_->data();
  for (auto& abit : det_->string_bits_a()) {
    for (auto& bbit : det_->string_bits_b()) {
      ordered_elements.emplace(*d++, make_pair(abit, bbit));
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
      const size_t ia = det_->lexical<0>(basis[i].first);
      const size_t ib = det_->lexical<1>(basis[i].second);

      for (int j = 0; j < nstate_; ++j)
        out->data(j)->element(ib, ia) = coeffs1->element(i, j);
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
void FCI::generate_guess(const int nspin, const int nstate, shared_ptr<Dvec> out) {
  int ndet = nstate_*10;
  start_over:
  vector<pair<bitset<nbit__>, bitset<nbit__>>> bits = detseeds(ndet);

  // Spin adapt detseeds
  int oindex = 0;
  vector<bitset<nbit__>> done_open;
  vector<bitset<nbit__>> done_closed;
  for (auto& it : bits) {
    bitset<nbit__> alpha = it.second;
    bitset<nbit__> beta = it.first;
    bitset<nbit__> open_bit = (alpha^beta);
    bitset<nbit__> closed_bit = (alpha&beta);

    // This can happen if all possible determinants are checked without finding nstate acceptable ones.
    if (alpha.count() + beta.count() != nelea_ + neleb_)
      throw logic_error("FCI::generate_guess produced an invalid determinant.  Check the number of states being requested.");

    // make sure that we have enough unpaired alpha
    const int unpairalpha = (alpha ^ (alpha & beta)).count();
    const int unpairbeta  = (beta ^ (alpha & beta)).count();
    if (unpairalpha-unpairbeta < nelea_-neleb_) continue;

    // check if this orbital configuration is already used
    if ((find(done_open.begin(), done_open.end(), open_bit) != done_open.end()) && (find(done_closed.begin(), done_closed.end(), closed_bit) != done_closed.end())) continue;
    done_open.push_back(open_bit);
    done_closed.push_back(closed_bit);

    pair<vector<tuple<int, int, int>>, double> adapt = det()->spin_adapt(nelea_-neleb_, alpha, beta);
    const double fac = adapt.second;
    for (auto& iter : adapt.first) {
      out->data(oindex)->element(get<0>(iter), get<1>(iter)) = get<2>(iter)*fac;
    }
    out->data(oindex)->spin_decontaminate();

    cout << "     guess " << setw(3) << oindex << ":   closed " <<
          setw(20) << left << print_bit(closed_bit, norb_) << " open " << setw(20) << print_bit(open_bit, norb_) << right << endl;

    ++oindex;
    if (oindex == nstate) break;
  }
  if (oindex < nstate) {
    out->zero();
    ndet *= 4;
    goto start_over;
  }
  cout << endl;
}

// returns seed determinants for initial guess
vector<pair<bitset<nbit__> , bitset<nbit__>>> FCI::detseeds(const int ndet) const {
  multimap<double, pair<bitset<nbit__>,bitset<nbit__>>> tmp;
  for (int i = 0; i != ndet; ++i) tmp.emplace(-1.0e10*(1+i), make_pair(bitset<nbit__>(0),bitset<nbit__>(0)));

  double* diter = denom_->data();
  for (auto& aiter : det()->string_bits_a()) {
    for (auto& biter : det()->string_bits_b()) {
      const double din = -(*diter);
      if (tmp.begin()->first < din) {
        tmp.emplace(din, make_pair(biter, aiter));
        tmp.erase(tmp.begin());
      }
      ++diter;
    }
  }
  assert(tmp.size() == ndet || ndet > det()->string_bits_a().size()*det()->string_bits_b().size());
  vector<pair<bitset<nbit__> , bitset<nbit__>>> out;
  for (auto iter = tmp.rbegin(); iter != tmp.rend(); ++iter)
    out.push_back(iter->second);
  return out;
}


void FCI::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "        FCI calculation      " << endl;
  cout << "  ---------------------------" << endl << endl;
}


shared_ptr<const CIWfn> FCI::conv_to_ciwfn() const {
  return make_shared<CIWfn>(geom_, ncore_, norb_, nstate_, energy_, cc_, det_);
}


void FCI::compute() {
  Timer pdebug(3);

  if (!restarted_) {
    // Creating an initial CI vector
    cc_ = make_shared<Dvec>(det_, nstate_); // B runs first

    // find determinants that have small diagonal energies
    if (nguess_ <= nstate_)
      generate_guess(nelea_-neleb_, nstate_, cc_);
    else
      model_guess(cc_);
    pdebug.tick_print("guess generation");

    // Davidson utility
    davidson_ = make_shared<DavidsonDiag<Civec>>(nstate_, davidson_subspace_);
  }

  // nuclear energy retrieved from geometry
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();

  // main iteration starts here
  cout << "  === FCI iteration ===" << endl << endl;
  // 0 means not converged
  vector<int> conv(nstate_, 0);

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer fcitime;

    // form a sigma vector given cc
    shared_ptr<Dvec> sigma = form_sigma(cc_, jop_, conv);
    pdebug.tick_print("sigma vector");

#ifndef DISABLE_SERIALIZATION
    if (restart_) {
      stringstream ss; ss << "fci_" << iter;
      OArchive ar(ss.str());
      ar << static_cast<Method*>(this);
    }
#endif

    // constructing Dvec's for Davidson
    auto ccn = make_shared<const CASDvec>(cc_->dvec());
    auto sigman = make_shared<const CASDvec>(sigma->dvec());
    const vector<double> energies = davidson_->compute(ccn->dvec(conv), sigman->dvec(conv));

    // get residual and new vectors
    vector<shared_ptr<Civec>> errvec = davidson_->residual();
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
        double* source_array = errvec[ist]->data();
        double* denom_array = denom_->data();
        const double en = energies[ist];
        for (size_t i = 0; i != size; ++i) {
          target_array[i] = source_array[i] / min(en - denom_array[i], -0.1);
        }
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

  auto cc = make_shared<CASDvec>(davidson_->civec());
  cc_ = make_shared<Dvec>(*cc);
  cc_->print(print_thresh_);

  if (dipoles_) {
    compute_rdm12();
    cout << endl;
    auto ocoeff = coeff_->slice(0, ncore_+norb_);
    for (int i = 0; i != nstate_; ++i) {
      stringstream ss; ss << "state " << i;
      shared_ptr<const Matrix> mordm = rdm1(i)->rdm1_mat(ncore_);
      Dipole dipole(geom_, make_shared<Matrix>(ocoeff * *mordm ^ ocoeff), ss.str());
      dipole.compute();
    }
    for (int i = 0; i != nstate_; ++i)
      for (int j = 0; j != i; ++j) {
        stringstream ss; ss << "states " << i << " " << j;
        compute_rdm12(i, j);
        shared_ptr<const Matrix> mordm = rdm1(i,j)->rdm1_mat(ncore_, false);
        Dipole dipole(geom_, make_shared<Matrix>(ocoeff * *mordm ^ ocoeff), ss.str());
        dipole.compute();
      }
  }
}
