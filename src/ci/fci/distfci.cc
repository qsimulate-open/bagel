//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci/distfci.cc
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

#include <cassert>

#include <src/util/combination.hpp>
#include <src/util/math/comb.h>
#include <src/util/math/davidson.h>
#include <src/ci/fci/modelci.h>
#include <src/ci/fci/dist_form_sigma.h>
#include <src/ci/fci/distfci.h>
#include <src/ci/fci/space.h>
#include <src/ci/fci/fci_base.h>
#include <src/ci/fci/hzdenomtask.h>

using namespace std;
using namespace bagel;


void DistFCI::common_init() {
#ifndef HAVE_MPI_H
  throw logic_error("DistFCI can be used only with MPI");
#endif

  cout << "    * Parallel algorithm will be used." << endl;

  print_header();

  const bool frozen = idata_->get<bool>("frozen", false);
  max_iter_ = idata_->get<int>("maxiter", 100);
  max_iter_ = idata_->get<int>("maxiter_fci", max_iter_);
  davidson_subspace_ = idata_->get<int>("davidson_subspace", 20);
  thresh_ = idata_->get<double>("thresh", 1.0e-10);
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
  weight_ = vector<double>(nstate_, 1.0/static_cast<double>(nstate_));

  rdm1_ = make_shared<VecRDM<1>>();
  rdm2_ = make_shared<VecRDM<2>>();
  energy_.resize(nstate_);

  // construct a determinant space in which this FCI will be performed.
  space_ = make_shared<HZSpace>(norb_, nelea_, neleb_);
  det_ = space_->finddet(nelea_, neleb_);
}


void DistFCI::model_guess(vector<shared_ptr<DistCivec>>& out) {
  multimap<double, pair<size_t, size_t>> ordered_elements;
  {
    const double* d = denom_->local_data();
    for (size_t ia = denom_->astart(); ia < denom_->aend(); ++ia)
      for (size_t ib = 0; ib < det_->lenb(); ++ib)
        ordered_elements.emplace(*d++, make_pair(ia, ib));
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
      basis.emplace_back(det_->string_bits_a(p.second.first), det_->string_bits_b(p.second.second));
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

    auto coeffs1 = (coeffs * *hamiltonian).slice_copy(0, nstate_);
    mpi__->broadcast(coeffs1->data(), coeffs1->ndim() * coeffs1->mdim(), 0);
    for (int i = 0; i < nguess; ++i) {
      size_t ia = det_->lexical<0>(basis[i].first);
      if ( ia >= denom_->astart() && ia < denom_->aend() ) {
        ia -= denom_->astart();
        const size_t ib = det_->lexical<1>(basis[i].second);
        for (int j = 0; j < nstate_; ++j)
          out[j]->set_local(ia, ib, coeffs1->element(i, j));
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
  vector<bitset<nbit__>> done_open;
  vector<bitset<nbit__>> done_closed;
  for (auto& it : bits) {
    bitset<nbit__> alpha = it.second;
    bitset<nbit__> beta = it.first;
    bitset<nbit__> open_bit = (alpha^beta);
    bitset<nbit__> closed_bit = (alpha&beta);

    // This can happen if all possible determinants are checked without finding nstate acceptable ones.
    if (alpha.count() + beta.count() != nelea_ + neleb_)
      throw logic_error("DistFCI::generate_guess produced an invalid determinant.  Check the number of states being requested.");

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

    out[oindex]->zero();
    for (auto& ad : adapt.first) {
      const int aloc = get<1>(ad) - out[oindex]->astart();
      if (aloc >= 0 && aloc < out[oindex]->asize())
        out[oindex]->set_local(aloc, get<0>(ad), get<2>(ad)*fac);
    }
    out[oindex]->spin_decontaminate();

    cout << "     guess " << setw(3) << oindex << ":   closed " <<
          setw(20) << left << print_bit(closed_bit, norb_) << " open " << setw(20) << print_bit(open_bit, norb_) << right << endl;

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

vector<pair<bitset<nbit__> , bitset<nbit__>>> DistFCI::detseeds(const int ndet) const {
  multimap<double, pair<size_t, size_t>> tmp;
  for (int i = 0; i != ndet; ++i)
    tmp.emplace(-1.0e10*(1+i), make_pair(0,0));

  const double* diter = denom_->local_data();
  for (size_t ia = denom_->astart(); ia != denom_->aend(); ++ia) {
    for (size_t ib = 0; ib != det_->lenb(); ++ib) {
      const double din = -*diter++;
      if (tmp.begin()->first < din) {
        tmp.emplace(din, make_pair(ib, ia));
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
    tmp.emplace(eall[i], make_pair(ball[i], aall[i]));
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
    out.push_back({det_->string_bits_b(ball[i]), det_->string_bits_a(aall[i])});

  return out;
}


void DistFCI::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "        FCI calculation      " << endl;
  cout << "  ---------------------------" << endl << endl;
}



void DistFCI::update(shared_ptr<const Matrix> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  coeff_ = make_shared<Matrix>(*c);
  jop_ = make_shared<Jop>(ref_, ncore_, ncore_+norb_, coeff_, store_half_ints_, "HZ");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;

  const_denom();
}



// same as HZ::const_denom except that denom_ is also distributed
void DistFCI::const_denom() {
  Timer denom_t;
  auto h = make_shared<VectorB>(norb_);
  auto jop = make_shared<Matrix>(norb_, norb_);
  auto kop = make_shared<Matrix>(norb_, norb_);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      (*jop)(i,j) = (*jop)(j,i) = 0.5*jop_->mo2e_hz(i, j, i, j);
      (*kop)(i,j) = (*kop)(j,i) = 0.5*jop_->mo2e_hz(i, j, j, i);
    }
    (*h)(i) = jop_->mo1e(i,i);
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<DistCivec>(det_);

  unique_ptr<double[]> buf(new double[denom_->size()]);
  auto iter = buf.get();
  TaskQueue<HZDenomTask> tasks(denom_->asize());
  for (size_t i = denom_->astart(); i != denom_->aend(); ++i) {
    tasks.emplace_back(iter, denom_->det()->string_bits_a(i), det_, jop, kop, h);
    iter += det()->string_bits_b().size();
  }
  tasks.compute();

  denom_->accumulate_buffer(1.0, buf);
  denom_t.tick_print("denom");
}


void DistFCI::compute() {
  Timer pdebug(3);

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
  DavidsonDiag<DistCivec> davidson(nstate_, davidson_subspace_);

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

    vector<shared_ptr<const DistCivec>> ccn, sigman;
    for (int i = 0; i < nstate_; ++i) {
      ccn.push_back(cc[i]);
      sigman.push_back(sigma[i]);
    }

    // Davidson
    const vector<double> energies = davidson.compute(ccn, sigman);

    // get residual and new vectors
    vector<shared_ptr<DistCivec>> errvec = davidson.residual();
    pdebug.tick_print("davidson");

    // compute errors
    vector<double> errors;
    for (int i = 0; i != nstate_; ++i) {
      errors.push_back(errvec[i]->rms());
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
          unique_ptr<double[]> target_array(new double[c->size()]);
          const double* source_array = errvec[ist]->local_data();
          const double* denom_array = denom_->local_data();
          const double en = energies[ist];
          // TODO this should be threaded
          for (int i = 0; i != size; ++i) {
            target_array[i] = source_array[i] / min(en - denom_array[i], -0.1);
          }
          c->accumulate_buffer(1.0, target_array);
          c->spin_decontaminate();
          c->normalize();
          cc.push_back(c);
        } else {
          cc.push_back(nullptr);
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

shared_ptr<const Civec> DistFCI::denom() const {
  return denom_->civec();
}

shared_ptr<const CIWfn> DistFCI::conv_to_ciwfn() const {
  auto cc = conv_to_dvec();
  return make_shared<CIWfn>(geom_, ncore_, norb_, nstate_, energy_, cc, cc->det());
}

shared_ptr<const Dvec> DistFCI::civectors() const {
  auto cc = conv_to_dvec();
  return cc;
}

shared_ptr<const Dvec> DistFCI::conv_to_dvec() const {
  vector<shared_ptr<Civec>> vec;
  for (auto& i : cc_->dvec())
    vec.push_back(i->civec());
  auto cc = make_shared<Dvec>(Dvector_base<Civec>(vec));
  return cc;
}

shared_ptr<Dvec> DistFCI::distdvec_to_dvec(shared_ptr<const DistDvec> d) const {
  vector<shared_ptr<Civec>> vec;
  for (auto& i : d->dvec())
    vec.push_back(i->civec());
  auto cc = make_shared<Dvec>(Dvector_base<Civec>(vec));
  return cc;
}

shared_ptr<const DistDvec> DistFCI::dvec_to_distdvec(shared_ptr<const Dvec> c) const {
  vector<shared_ptr<DistCivec>> vec;
  for (auto& i : c->dvec())
    vec.push_back(make_shared<DistCivec>(i));
  auto cc = make_shared<DistDvec>(Dvector_base<DistCivec>(vec));
  return cc;
}

tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
DistFCI::compute_rdm12_av_from_dvec(shared_ptr<const Dvec> dbrab, shared_ptr<const Dvec> dketb, shared_ptr<const Determinants> o) const {
  auto dbra = dvec_to_distdvec(dbrab);
  auto dket = dvec_to_distdvec(dketb);

  if (o != nullptr) {
    dbra->set_det(o);
    dket->set_det(o);
  }

  auto rdm1 = make_shared<RDM<1>>(norb_);
  auto rdm2 = make_shared<RDM<2>>(norb_);

  assert(dbra->ij() == dket->ij() && dbra->det() == dket->det());

  for (int i = 0; i != dbra->ij(); ++i) {
    shared_ptr<RDM<1>> r1;
    shared_ptr<RDM<2>> r2;
    tie(r1, r2) = compute_rdm12_from_civec(dbra->data(i), dket->data(i));
    rdm1->ax_plus_y(weight_[i], r1);
    rdm2->ax_plus_y(weight_[i], r2);
  }

  if (o != nullptr) {
    dbra->set_det(det_);
    dket->set_det(det_);
  }

  return tie(rdm1, rdm2);
}
