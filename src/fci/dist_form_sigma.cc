//
// BAGEL - Parallel electron correlation program.
// Filename: fci/dist_form_sigma.cc
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

#include <src/util/combination.hpp>
#include <src/math/comb.h>
#include <src/fci/distfci_ab.h>
#include <src/fci/distfci_bb.h>
#include <src/fci/dist_form_sigma.h>

using namespace std;
using namespace bagel;

vector<shared_ptr<DistCivec>> FormSigmaDistFCI::operator()(const vector<shared_ptr<DistCivec>>& ccvec, shared_ptr<const MOFile> jop, const vector<int>& conv) const {
  const int norb = jop->nocc();
  const int ij = norb*norb;

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
    sigma_aa(ctrans, strans, jop, cc->det()->remalpha()->rembeta());
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

shared_ptr<DistDvec> FormSigmaDistFCI::operator()(shared_ptr<const DistDvec> ccvec, shared_ptr<const MOFile> jop) const {
  vector<int> conv(ccvec->ij(), static_cast<bool>(false));
  vector<shared_ptr<DistCivec>> svec = (*this)(ccvec->dvec(), jop, conv);
  return make_shared<DistDvec>(svec);
}


void FormSigmaDistFCI::sigma_ab(shared_ptr<const DistCivec> cc, shared_ptr<DistCivec> sigma, shared_ptr<const MOFile> jop) const {
  shared_ptr<const Determinants> base_det = cc->det();
  shared_ptr<const Determinants> int_det = base_det->remalpha()->rembeta();

  const int norb = base_det->norb();

  const size_t lbt = int_det->lenb();
  const size_t lbs = base_det->lenb();
  const int ij = norb*norb;

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



void FormSigmaDistFCI::sigma_aa(shared_ptr<const DistCivec> ctrans, shared_ptr<DistCivec> strans, shared_ptr<const MOFile> jop, shared_ptr<const Determinants> int_det) const {
  shared_ptr<const Determinants> trans_det = ctrans->det();
  shared_ptr<const Determinants> int_tra = int_det->transpose();
  sigma_bb(ctrans, strans, jop, ctrans->det(), int_tra);
}


void FormSigmaDistFCI::sigma_bb(shared_ptr<const DistCivec> cc, shared_ptr<DistCivec> sigma, shared_ptr<const MOFile> jop) const {
  const shared_ptr<const Determinants> base_det = cc->det();
  const shared_ptr<const Determinants> int_det = base_det->remalpha()->rembeta(); // only for n-1 beta strings...
  sigma_bb(cc, sigma, jop, base_det, int_det);
}


// beta-beta block has no communication (and should be cheap)
void FormSigmaDistFCI::sigma_bb(shared_ptr<const DistCivec> cc, shared_ptr<DistCivec> sigma, shared_ptr<const MOFile> jop,
                       const shared_ptr<const Determinants> base_det, const shared_ptr<const Determinants> int_det) const {

  const int norb = cc->det()->norb();
  const size_t lb = sigma->lenb();
  const size_t la = sigma->asize();

  unique_ptr<double[]> target(new double[la*lb]);
  unique_ptr<double[]> source(new double[la*lb]);

  // (astart:aend, b)
  blas::transpose(cc->local(), lb, la, source.get());
  fill_n(target.get(), la*lb, 0.0);

  // preparing Hamiltonian
  const size_t npack = norb*(norb-1)/2;
  unique_ptr<double[]> hamil1(new double[norb*norb]);
  unique_ptr<double[]> hamil2(new double[npack*npack]);
  for (int i = 0, ij = 0, ijkl = 0; i != norb; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      hamil1[j+norb*i] = hamil1[i+norb*j] = jop->mo1e(ij);
      if (i == j) continue;
      for (int k = 0; k != norb; ++k)
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
  const size_t lengb = comb.c(norb, neleb-2);
  vector<bitset<nbit__>> intb(lengb, bitset<nbit__>(0));
  vector<int> data(norb);
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

  blas::transpose(target.get(), la, lb, source.get());
  for (size_t i = 0; i != la; ++i) {
    lock_guard<mutex> lock(sigma->cimutex(i));
    daxpy_(lb, 1.0, source.get()+i*lb, 1, sigma->local()+i*lb, 1);
  }

  // for accumulate in aa and ab
#ifndef USE_SERVER_THREAD
  sigma->flush();
#endif
}
