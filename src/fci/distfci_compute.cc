//
// BAGEL - Parallel electron correlation program.
// Filename: distfci.cc
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

#include <src/fci/distfci.h>
#include <src/math/davidson.h>
#include <src/util/combination.hpp>
#include <src/math/comb.h>
#include <src/fci/distfci_ab.h>
#include <src/fci/distfci_bb.h>

using namespace std;
using namespace bagel;

vector<shared_ptr<DistCivec>> DistFCI::form_sigma(vector<shared_ptr<DistCivec>>& ccvec, shared_ptr<const MOFile> jop, const vector<int>& conv) const {
  const int ij = norb_*norb_;

  const int nstate = ccvec.size();

  vector<shared_ptr<DistCivec>> sigmavec;

  for (int istate = 0; istate != nstate; ++istate) {
    if (conv[istate]) {
      sigmavec.push_back(shared_ptr<DistCivec>());
      continue;
    }
    shared_ptr<const DistCivec> cc = ccvec[istate];
    shared_ptr<DistCivec> sigma = cc->clone();
    sigma->zero();

    Timer fcitime(1);
    sigma->init_mpi_accumulate();

    shared_ptr<DistCivec> ctrans = cc->transpose();

    sigma_ab(cc, sigma, jop);
    fcitime.tick_print("alpha-beta");

    ctrans->transpose_wait();
    shared_ptr<DistCivec> strans = ctrans->clone();
    sigma_aa(ctrans, strans, jop);
    shared_ptr<DistCivec> sigma_aa = strans->transpose();
    fcitime.tick_print("alpha-alpha");

    sigma_bb(cc, sigma, jop);
    fcitime.tick_print("beta-beta");

    sigma_aa->transpose_wait();
    sigma->ax_plus_y(1.0, *sigma_aa);
    fcitime.tick_print("wait1");

    sigma->terminate_mpi_accumulate();
    fcitime.tick_print("wait");

    sigmavec.push_back(sigma);
  }

  return sigmavec;
}


void DistFCI::sigma_ab(shared_ptr<const DistCivec> cc, shared_ptr<DistCivec> sigma, shared_ptr<const MOFile> jop) const {
  shared_ptr<const Determinants> base_det = cc->det();
  shared_ptr<const Determinants> int_det = space_->finddet(base_det->nelea() - 1, base_det->neleb() - 1);

  const size_t lbt = int_det->lenb();
  const size_t lbs = base_det->lenb();
  const int ij = norb_*norb_;

  const int rank = mpi__->rank();
  const int size = mpi__->size();

  cc->init_mpi_recv();

  vector<shared_ptr<DistABTask>> tasks;

  // shamelessly statically distributing across processes
  for (size_t a = 0; a != int_det->lena(); ++a) {
    if (a%size != rank) continue;

    const bitset<nbit__> astring = int_det->stringa(a);

    tasks.push_back(make_shared<DistABTask>(astring, base_det, int_det, jop, cc, sigma));

    for (auto i = tasks.begin(); i != tasks.end(); ) {
      if ((*i)->test()) {
         (*i)->compute();
        i = tasks.erase(i);
      } else {
        ++i;
      }
    }
#ifndef USE_SERVER_THREAD
    cc->flush();
    sigma->flush();
#endif
  }

  bool done;
  do {
    done = true;
    for (auto i = tasks.begin(); i != tasks.end(); ) {
      if ((*i)->test()) {
        (*i)->compute();
        i = tasks.erase(i);
      } else {
        ++i;
        done = false;
      }
    }
#ifndef USE_SERVER_THREAD
    size_t d = done ? 0 : 1;
    mpi__->soft_allreduce(&d, 1);
    done = d == 0;
    if (!done) cc->flush();
    if (!done) sigma->flush();
#endif
    if (!done) this_thread::sleep_for(sleeptime__);
  } while (!done);

  cc->terminate_mpi_recv();
}



void DistFCI::sigma_aa(shared_ptr<const DistCivec> ctrans, shared_ptr<DistCivec> strans, shared_ptr<const MOFile> jop) const {
  shared_ptr<const Determinants> trans_det = ctrans->det();
  shared_ptr<const Determinants> int_tra = space_->finddet(trans_det->neleb() - 1, trans_det->nelea() - 1)->transpose();
  sigma_bb(ctrans, strans, jop, ctrans->det(), int_tra);
}


void DistFCI::sigma_bb(shared_ptr<const DistCivec> cc, shared_ptr<DistCivec> sigma, shared_ptr<const MOFile> jop) const {
  const shared_ptr<const Determinants> base_det = cc->det();
  const shared_ptr<const Determinants> int_det = space_->finddet(base_det->nelea() - 1,base_det->neleb() - 1); // only for n-1 beta strings...
  sigma_bb(cc, sigma, jop, base_det, int_det);
}


// beta-beta block has no communication (and should be cheap)
void DistFCI::sigma_bb(shared_ptr<const DistCivec> cc, shared_ptr<DistCivec> sigma, shared_ptr<const MOFile> jop,
                       const shared_ptr<const Determinants> base_det, const shared_ptr<const Determinants> int_det) const {

  const size_t lb = sigma->lenb();
  const size_t la = sigma->asize();

  unique_ptr<double[]> target(new double[la*lb]);
  unique_ptr<double[]> source(new double[la*lb]);

  // (astart:aend, b)
  mytranspose_(cc->local(), lb, la, source.get());
  fill_n(target.get(), la*lb, 0.0);

  // preparing Hamiltonian
  const size_t npack = norb_*(norb_-1)/2;
  unique_ptr<double[]> hamil1(new double[norb_*norb_]);
  unique_ptr<double[]> hamil2(new double[npack*npack]);
  for (int i = 0, ij = 0, ijkl = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      hamil1[j+norb_*i] = hamil1[i+norb_*j] = jop_->mo1e(ij);
      if (i == j) continue;
      for (int k = 0; k != norb_; ++k)
        for (int l = 0; l < k; ++l, ++ijkl)
          hamil2[ijkl] = jop->mo2e_hz(l,k,j,i) - jop->mo2e_hz(k,l,j,i);
    }
  }

#ifndef USE_SERVER_THREAD
  // for accumulate in aa and ab
  sigma->flush();
#endif

  const size_t nelea = base_det->nelea();
  const size_t neleb = base_det->neleb();

  const static Comb comb;
  const size_t lengb = comb.c(norb_, neleb-2);
  vector<bitset<nbit__>> intb(lengb, bitset<nbit__>(0));
  vector<int> data(norb_);
  iota(data.begin(), data.end(), 0);
  auto sa = intb.begin();
  do {
    for (int i=0; i < neleb-2; ++i) sa->set(data[i]);
    ++sa;
  } while (boost::next_combination(data.begin(), data.begin()+neleb-2, data.end()));

  vector<mutex> localmutex(lb);
  // loop over intermediate string
  TaskQueue<DistBBTask> tasks(intb.size());

  // two electron part
  for (auto& b : intb)
    tasks.emplace_back(la, source.get(), target.get(), hamil2.get(), base_det, b, &localmutex);
  // one electron part
  for (auto& b : int_det->stringb())
    tasks.emplace_back(la, source.get(), target.get(), hamil1.get(), base_det, b, &localmutex);

  tasks.compute();

  mytranspose_(target.get(), la, lb, source.get());
  for (size_t i = 0; i != la; ++i) {
    lock_guard<mutex> lock(sigma->cimutex(i));
    daxpy_(lb, 1.0, source.get()+i*lb, 1, sigma->local()+i*lb, 1);
  }

  // for accumulate in aa and ab
#ifndef USE_SERVER_THREAD
  sigma->flush();
#endif
}



void DistFCI::compute() {
  Timer pdebug(2);

  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("FCI: C1 only at the moment.");

  // some constants
  //const int ij = nij();

  // Creating an initial CI vector
  vector<shared_ptr<DistCivec>> cc(nstate_);
  for (auto& i : cc)
    i = make_shared<DistCivec>(det_);

  // find determinants that have small diagonal energies
  generate_guess(nelea_-neleb_, nstate_, cc);
  pdebug.tick_print("guess generation");

  // nuclear energy retrieved from geometry
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();

  // Davidson utility
  DavidsonDiag<DistCivec> davidson(nstate_, max_iter_);

  // main iteration starts here
  cout << "  === FCI iteration ===" << endl << endl;
  // 0 means not converged
  vector<int> conv(nstate_,0);

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
