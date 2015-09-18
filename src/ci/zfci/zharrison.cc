//
// BAGEL - Parallel electron correlation program.
// Filename: zharrison.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/ci/zfci/zharrison.h>
#include <src/ci/zfci/relspace.h>
#include <src/util/math/comb.h>
#include <src/prop/pseudospin/pseudospin.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/mat1e/giao/relhcore_london.h>
#include <src/mat1e/rel/reloverlap.h>
#include <src/mat1e/giao/reloverlap_london.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::ZHarrison)

using namespace std;
using namespace bagel;

ZHarrison::ZHarrison(std::shared_ptr<const PTree> idat, shared_ptr<const Geometry> g, shared_ptr<const Reference> r, const int ncore, const int norb, const int nstate, std::shared_ptr<const RelCoeff_Block> coeff_zcas, const bool restricted)
 : Method(idat, g, r), ncore_(ncore), norb_(norb), nstate_(nstate), restarted_(false) {
  if (!ref_) throw runtime_error("ZFCI requires a reference object");

  auto rr = dynamic_pointer_cast<const RelReference>(ref_);
  if (!rr) throw runtime_error("ZFCI currently requires a relativistic reference object");

  const bool frozen = idata_->get<bool>("frozen", false);
  max_iter_ = idata_->get<int>("maxiter", 100);
  max_iter_ = idata_->get<int>("maxiter_fci", max_iter_);
  davidson_subspace_ = idata_->get<int>("davidson_subspace", 20);
  thresh_ = idata_->get<double>("thresh", 1.0e-10);
  thresh_ = idata_->get<double>("thresh_fci", thresh_);
  print_thresh_ = idata_->get<double>("print_thresh", 0.05);
  restart_ = idata_->get<bool>("restart", false);

  states_ = idata_->get_vector<int>("state", 0);
  nstate_ = 0;
  for (int i = 0; i != states_.size(); ++i)
    nstate_ += states_[i] * (i+1); // 2S+1 for 0, 1/2, 1, ...

  // if nstate was specified on construction, it should match
  assert(nstate == nstate_ || nstate == -1);

  gaunt_ = idata_->get<bool>("gaunt", rr->gaunt());
  breit_ = idata_->get<bool>("breit", rr->breit());
  if (gaunt_ != rr->gaunt())
    geom_ = geom_->relativistic(gaunt_);

  // Invoke Kramer's symmetry for any case without magnetic field
  tsymm_ = !geom_->magnetism();

  if (ncore_ < 0)
    ncore_ = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (norb_  < 0)
    norb_ = idata_->get<int>("norb", (rr->relcoeff()->mdim()/2-ncore_));

  // additional charge
  charge_ = idata_->get<int>("charge", 0);

  nele_ = geom_->nele() - charge_ - ncore_*2;

  if (norb_ < 0 || norb_ + ncore_ > geom_->nbasis())
    throw runtime_error("Invalid number of active orbitals");
  if (nele_ < 0)
    throw runtime_error("Number of active electrons is less than zero.");

  energy_.resize(nstate_);

  print_header();
  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << ncore_ << endl;
  cout << "    * nact     : " << setw(6) << norb_ << endl;
  cout << "    * nvirt    : " << setw(6) << ((coeff_zcas ? coeff_zcas->mdim() : rr->relcoeff_full()->mdim())/2-ncore_-norb_) << endl;

  if (!geom_->dfs())
    geom_ = geom_->relativistic(gaunt_);

  space_ = make_shared<RelSpace>(norb_, nele_);
  int_space_ = make_shared<RelSpace>(norb_, nele_-2, /*mute*/true, /*link up*/true);
  assert((restricted && coeff_zcas) || (!restricted && !coeff_zcas));

  // obtain the coefficient matrix in striped format
  shared_ptr<const RelCoeff_Block> coeff;
  if (coeff_zcas) {
    coeff = coeff_zcas;
  } else {
    auto scoeff = make_shared<const RelCoeff_Striped>(*rr->relcoeff_full(), ncore_, norb_, rr->relcoeff_full()->mdim()/4-ncore_-norb_, rr->relcoeff_full()->mdim()/2);

    shared_ptr<const ZMatrix> hcore, overlap;
    if (!geom_->magnetism()) {
      overlap = make_shared<RelOverlap>(geom_);
      hcore = make_shared<RelHcore>(geom_);
    } else {
      overlap = make_shared<RelOverlap_London>(geom_);
      hcore = make_shared<RelHcore_London>(geom_);
    }

    // then compute Kramers adapated coefficient matrices
    scoeff = scoeff->init_kramers_coeff(geom_, overlap, hcore, geom_->nele()-charge_, tsymm_, gaunt_, breit_);

    // generate modified virtual orbitals, if requested
    const bool mvo = idata_->get<bool>("generate_mvo", false);
    if (mvo) {
      const bool hcore_mvo = idata_->get<bool>("hcore_mvo", false);
      const int ncore_mvo = idata_->get<int>("ncore_mvo", geom_->num_count_ncore_only());
      if (ncore_mvo == 2*rr->relcoeff_full()->nocc())
        cout << "    +++ Modified virtuals are Dirac-Fock orbitals with this choice of the core +++ "<< endl;
      else
        scoeff = scoeff->generate_mvo(geom_, overlap, hcore, ncore_mvo, rr->relcoeff_full()->nocc(), hcore_mvo, tsymm_, gaunt_, breit_);
    }

    // Reorder as specified in the input so frontier orbitals contain the desired active space
    const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
    if (iactive) {
      set<int> active_indices;
      // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
      for (auto& i : *iactive)
        active_indices.insert(lexical_cast<int>(i->data()) - 1);
      scoeff = scoeff->set_active(active_indices, geom_->nele()-charge_, tsymm_);
    }
    coeff = scoeff->block_format();
  }

  update(coeff, restricted);

}


void ZHarrison::print_header() const {
  cout << "  ----------------------------" << endl;
  cout << "  Relativistic FCI calculation" << endl;
  cout << "  ----------------------------" << endl << endl;
  cout << "    * Correlation of " << nele_ << " active electrons in " << norb_ << " orbitals."  << endl;
  cout << "    * Time-reversal symmetry " << (tsymm_ ? "will be assumed." : "violation will be permitted.") << endl;
  cout << "    * gaunt    : " << (gaunt_ ? "true" : "false") << endl;
  cout << "    * breit    : " << (breit_ ? "true" : "false") << endl;
}


// generate initial vectors
void ZHarrison::generate_guess(const int nelea, const int neleb, const int nstate, std::shared_ptr<RelZDvec> out, const int offset) {
  shared_ptr<const Determinants> cdet = space_->finddet(nelea, neleb);
  int ndet = nstate*10;
  int oindex = offset;
  const bool spin_adapt = idata_->get<bool>("spin_adapt", true);
  while (oindex < offset+nstate) {
    vector<pair<bitset<nbit__>, bitset<nbit__>>> bits = detseeds(ndet, nelea, neleb);

    // Spin adapt detseeds
    oindex = offset;
    vector<pair<bitset<nbit__>,bitset<nbit__>>> done;
    for (auto& it : bits) {
      bitset<nbit__> alpha = it.second;
      bitset<nbit__> beta = it.first;
      bitset<nbit__> open_bit = (alpha^beta);

      // This can happen if all possible determinants are checked without finding nstate acceptable ones.
      if (alpha.count() + beta.count() != nele_)
        throw logic_error("ZFCI::generate_guess produced an invalid determinant.  Check the number of states being requested.");

      pair<bitset<nbit__>,bitset<nbit__>> config = spin_adapt ? make_pair(open_bit, alpha & beta) : it;
      if (find(done.begin(), done.end(), config) != done.end()) continue;
      done.push_back(config);

      // make sure that we have enough unpaired alpha
      const int unpairalpha = (alpha ^ (alpha & beta)).count();
      const int unpairbeta  = (beta ^ (alpha & beta)).count();
      if (unpairalpha-unpairbeta < nelea-neleb) continue;

      //if (find(done.begin(), done.end(), open_bit) != done.end()) continue;

      //done.push_back(open_bit);
      pair<vector<tuple<int, int, int>>, double> adapt;
       if (spin_adapt) {
         adapt = space_->finddet(nelea, neleb)->spin_adapt(nelea-neleb, alpha, beta);
       } else {
         adapt.first = vector<tuple<int, int, int>>(1, make_tuple(space_->finddet(nelea, neleb)->lexical<1>(beta),
                                                                  space_->finddet(nelea, neleb)->lexical<0>(alpha), 1));
         adapt.second = 1.0;
       }

      const double fac = adapt.second;
      for (auto& iter : adapt.first) {
        out->find(nelea, neleb)->data(oindex)->element(get<0>(iter), get<1>(iter)) = get<2>(iter)*fac;
      }
      cout << "     guess " << setw(3) << oindex << ":   closed " <<
            setw(20) << left << print_bit(alpha&beta, norb_) << " open " << setw(20) << print_bit(open_bit, norb_) << right << endl;

      ++oindex;
      if (oindex == offset+nstate) break;
    }

    if (oindex < offset+nstate) {
      for (int i = offset; i != offset+oindex; ++i) {
        out->find(nelea, neleb)->data(i)->zero();
      }
      ndet *= 4;
    }
  }
  assert(oindex == offset+nstate);
  cout << endl;
}


// returns seed determinants for initial guess
vector<pair<bitset<nbit__> , bitset<nbit__>>> ZHarrison::detseeds(const int ndet, const int nelea, const int neleb) const {
  shared_ptr<const Determinants> cdet = space_->finddet(nelea, neleb);

  multimap<double, pair<bitset<nbit__>,bitset<nbit__>>> tmp;
  for (int i = 0; i != ndet; ++i) tmp.emplace(-1.0e10*(1+i), make_pair(bitset<nbit__>(0),bitset<nbit__>(0)));

  double* diter = denom_->find(cdet->nelea(), cdet->neleb())->data();
  for (auto& aiter : cdet->string_bits_a()) {
    for (auto& biter : cdet->string_bits_b()) {
      const double din = -(*diter);
      if (tmp.begin()->first < din) {
        tmp.emplace(din, make_pair(biter, aiter));
        tmp.erase(tmp.begin());
      }
      ++diter;
    }
  }
  assert(tmp.size() == ndet || ndet > cdet->string_bits_a().size()*cdet->string_bits_b().size());
  vector<pair<bitset<nbit__> , bitset<nbit__>>> out;
  for (auto iter = tmp.rbegin(); iter != tmp.rend(); ++iter)
    out.push_back(iter->second);
  return out;
}


void ZHarrison::compute() {
  Timer pdebug(2);

  if (geom_->nirrep() > 1) throw runtime_error("ZFCI: C1 only at the moment.");

  if (!restarted_) {
    // Creating an initial CI vector
    cc_ = make_shared<RelZDvec>(space_, nstate_); // B runs first

    // TODO really we should check the number of states for each S value, rather than total number
    const static Comb combination;
    const size_t max_states = combination(2*norb_, nele_);
    if (nstate_ > max_states) {
      const string space = "(" + to_string(nele_) + "," + to_string(norb_) + ")";
      throw runtime_error("Wrong states specified - a " + space + " active space can only produce " + to_string(max_states) + " eigenstates.");
    }

    // find determinants that have small diagonal energies
    int offset = 0;
    for (int ispin = 0; ispin != states_.size(); ++ispin) {
      int nstate = 0;
      for (int i = ispin; i != states_.size(); ++i)
        nstate += states_[i];

      if (nstate == 0)
        continue;

      if ((geom_->nele()+ispin-charge_) % 2 == 1) {
        if (states_[ispin] == 0) {
          continue;
        } else {
          if ((geom_->nele()-charge_) % 2 == 0) throw runtime_error("Wrong states specified - only integer spins are allowed for even electron counts.");
          else throw runtime_error("Wrong states specified - only half-integer spins are allowed for odd electron counts.");
        }
      }

      const int nelea = (geom_->nele()+ispin-charge_)/2 - ncore_;
      const int neleb = (geom_->nele()-ispin-charge_)/2 - ncore_;
      if (neleb < 0) throw runtime_error("Wrong states specified - there are not enough active electrons for the requested spin state.");
      if (nelea > norb_) throw runtime_error("Wrong states specified - there are not enough active orbitals for the requested spin state.");

      generate_guess(nelea, neleb, nstate, cc_, offset);
      offset += nstate;
      if (nelea != neleb) {
        generate_guess(neleb, nelea, nstate, cc_, offset);
        offset += nstate;
      }
    }
    pdebug.tick_print("guess generation");

    // Davidson utility
    davidson_ = make_shared<DavidsonDiag<RelZDvec, ZMatrix>>(nstate_, davidson_subspace_);
  }

  // nuclear energy retrieved from geometry
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();

  // main iteration starts here
  cout << "  === Relativistic FCI iteration ===" << endl << endl;
  // 0 means not converged
  vector<int> conv(nstate_,0);

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer fcitime;

#ifndef DISABLE_SERIALIZATION
    if (restart_) {
      stringstream ss; ss << "zfci_" << iter;
      OArchive ar(ss.str());
      ar << static_cast<Method*>(this);
    }
#endif

    // form a sigma vector given cc
    shared_ptr<RelZDvec> sigma = form_sigma(cc_, jop_, conv);
    pdebug.tick_print("sigma vector");

    const vector<double> energies = davidson_->compute(cc_->dvec(conv), sigma->dvec(conv));
    // get residual and new vectors
    vector<shared_ptr<RelZDvec>> errvec = davidson_->residual();
    for (auto& i : errvec)
      i->synchronize();
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

      auto ctmp = errvec.front()->clone();

      for (int ist = 0; ist != nstate_; ++ist) {
        if (conv[ist]) continue;
        for (auto& ib : space_->detmap()) {
          const int na = ib.second->nelea();
          const int nb = ib.second->neleb();
          const size_t size = cc_->find(na, nb)->data(ist)->size();
          complex<double>* target_array = ctmp->find(na, nb)->data();
          complex<double>* source_array = errvec[ist]->find(na, nb)->data();
          double* denom_array = denom_->find(na, nb)->data();
          const double en = energies[ist];
          for (int i = 0; i != size; ++i) {
            target_array[i] = source_array[i] / min(en - denom_array[i], -0.1);
          }
        }
        ctmp->normalize();
        cc_->set_data(ist, ctmp);
      }
    }
    pdebug.tick_print("denominator");

    // printing out
    if (nstate_ != 1 && iter) cout << endl;
    for (int i = 0; i != nstate_; ++i) {
      cout << setw(7) << iter << setw(4) << i << " " << setw(2) << (conv[i] ? "*" : " ")
                              << setw(17) << fixed << setprecision(8) << energies[i]+nuc_core << "   "
                              << setw(10) << scientific << setprecision(2) << errors[i] << fixed << setw(10) << setprecision(2)
                              << fcitime.tick() << endl;
      energy_[i] = energies[i]+nuc_core;
    }
    if (*min_element(conv.begin(), conv.end())) break;
  }
  // main iteration ends here

  cc_ = make_shared<RelZDvec>(davidson_->civec());
  cc_->print(print_thresh_);

#if 0
  for (auto& iprop : properties_) {
    iprop->compute(cc_);
    iprop->print();
  }
#endif


  // TODO When the Property class is implemented, this should be one
  if (idata_->get<bool>("aniso", false)) {
    if (geom_->magnetism()) {
      cout << "  ** Magnetic anisotropy analysis is currently only available for zero-field calculations; sorry." << endl;
    } else {
      int nspin = idata_->get<int>("aniso_spin", states_.size()-1);
      Pseudospin ps(nspin);

      ps.compute_numerical_hamiltonian(*this, jop_->coeff_input()->active_part());

      vector<Spin_Operator> Dop = ps.build_2ndorder_zfs_operators();

      const bool real = idata_->get<bool>("aniso_real", false);
      Dop = ps.extract_hamiltonian_parameters(real, Dop);
      shared_ptr<const ZMatrix> D_1 = ps.compute_Dtensor(Dop);

      { // Extended Stevens Operators: default should grab the nonzero time-reversal symmetric orders, but can be specified in input
        vector<int> ranks = {};
        for (int i = 2; i <= nspin; i += 2)
          ranks.push_back(i);

        const shared_ptr<const PTree> eso_ranks = idata_->get_child_optional("aniso_rank");
        if (eso_ranks) {
          ranks = {};
          for (auto& i : *eso_ranks)
            ranks.push_back(lexical_cast<int>(i->data()));
        }
        vector<Spin_Operator> ESO = ps.build_extended_stevens_operators(ranks);
        ESO = ps.extract_hamiltonian_parameters(real, ESO);
        shared_ptr<const ZMatrix> D_2 = ps.compute_Dtensor(Dop);
        D_1->print("D-tensor, obtained directly");
        D_2->print("D-tensor, obtained from ESO");
      }
    }
  }
}


shared_ptr<const ZMatrix> ZHarrison::swap_pos_neg(shared_ptr<const ZMatrix> coeffin) const {
  auto out = coeffin->clone();
  const int n = coeffin->ndim();
  const int m = coeffin->mdim()/2;
  assert(n % 4 == 0 && m % 2 == 0 && m * 2 == coeffin->mdim());
  out->copy_block(0, 0, n, m, coeffin->get_submatrix(0, m, n, m));
  out->copy_block(0, m, n, m, coeffin->get_submatrix(0, 0, n, m));
  return out;
}


shared_ptr<const RelCIWfn> ZHarrison::conv_to_ciwfn() const {
  using PairType = pair<shared_ptr<const RelSpace>,shared_ptr<const RelSpace>>;
  return make_shared<RelCIWfn>(geom_, ncore_, norb_, nstate_, energy_, cc_, make_shared<PairType>(make_pair(space_, int_space_)));
}


