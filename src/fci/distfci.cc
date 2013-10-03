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

#include <src/fci/distfci.h>
#include <src/fci/space.h>
#include <src/util/combination.hpp>
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

  space_ = make_shared<Space>(det_, 1); 
  update(ref_->coeff());
}

void DistFCI::common_init() {
  print_header();

  const bool frozen = idata_->get<bool>("frozen", false);
  max_iter_ = idata_->get<int>("maxiter", 100);
  max_iter_ = idata_->get<int>("maxiter_fci", max_iter_);
  thresh_ = idata_->get<double>("thresh", 1.0e-20);
  thresh_ = idata_->get<double>("thresh_fci", thresh_);
  print_thresh_ = idata_->get<double>("print_thresh", 0.05);

  if (nstate_ < 0) nstate_ = idata_->get<int>("nstate", 1);

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
  det_ = make_shared<const Determinants>(norb_, nelea_, neleb_);
}

// generate initial vectors
//   - bits: bit patterns of low-energy determinants
//   - nspin: #alpha - #beta
//   - out:
void DistFCI::generate_guess(const int nspin, const int nstate, vector<shared_ptr<DistCivec>> out) {
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
  mpi__->allgather(&aarray[0], ndet, &aall[0], ndet);
  mpi__->allgather(&barray[0], ndet, &ball[0], ndet);
  mpi__->allgather(&en[0],     ndet, &eall[0], ndet);

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
  mpi__->broadcast(&aall[0], ndet, 0);
  mpi__->broadcast(&ball[0], ndet, 0);

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
  unique_ptr<double[]> h(new double[norb_]);
  unique_ptr<double[]> jop(new double[norb_*norb_]);
  unique_ptr<double[]> kop(new double[norb_*norb_]);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*jop_->mo2e_hz(i, j, i, j);
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*jop_->mo2e_hz(i, j, j, i);
    }
    h[i] = jop_->mo1e(i,i);
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<DistCivec>(det_);

  double* iter = denom_->local();
  TaskQueue<HZDenomTask> tasks(denom_->asize());
  for (size_t i = denom_->astart(); i != denom_->aend(); ++i) {
    tasks.emplace_back(iter, denom_->det()->stringa(i), det_, jop.get(), kop.get(), h.get());
    iter += det()->stringb().size();
  }
  tasks.compute();

  denom_t.tick_print("denom");
}
