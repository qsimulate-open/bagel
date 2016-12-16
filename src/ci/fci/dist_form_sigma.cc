//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci/dist_form_sigma.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/util/combination.hpp>
#include <src/util/math/comb.h>
#include <src/ci/fci/distfci_ab.h>
#include <src/ci/fci/distfci_bb.h>
#include <src/ci/fci/dist_form_sigma.h>

using namespace std;
using namespace bagel;

vector<shared_ptr<DistCivec>> FormSigmaDistFCI::operator()(const vector<shared_ptr<DistCivec>>& ccvec, shared_ptr<const MOFile> jop, const vector<int>& conv) const {
  const int nstate = ccvec.size();

  vector<shared_ptr<DistCivec>> sigmavec;

  for (int istate = 0; istate != nstate; ++istate) {
    if (conv[istate]) {
      sigmavec.push_back(nullptr);
      continue;
    }
    shared_ptr<const DistCivec> cc = ccvec[istate];
    shared_ptr<DistCivec> sigma = cc->clone();
    sigma->zero();

    Timer fcitime(3);

    sigma_ab(cc, sigma, jop);
    fcitime.tick_print("alpha-beta");

    shared_ptr<DistCivec> ctrans = cc->transpose();
    shared_ptr<DistCivec> strans = ctrans->clone();
    sigma_aa(ctrans, strans, jop, cc->det()->remalpha()->rembeta());
    fcitime.tick_print("alpha-alpha");

    sigma_bb(cc, sigma, jop);
    fcitime.tick_print("beta-beta");

    shared_ptr<DistCivec> sigma_aa = strans->transpose();
    sigma->ax_plus_y(1.0, *sigma_aa);
    fcitime.tick_print("post process");

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

  const int rank = mpi__->rank();
  const int size = mpi__->size();

  list<shared_ptr<DistABTask>> tasks;

  // shamelessly statically distributing across processes
  for (size_t a = 0; a != int_det->lena(); ++a) {
    if (a%size != rank) continue;

    const bitset<nbit__> astring = int_det->string_bits_a(a);
    tasks.push_back(make_shared<DistABTask>(astring, base_det, int_det, jop, cc, sigma));
  }

  list<shared_ptr<RMATask<double>>> acctasks;
  for (auto i = tasks.begin(); i != tasks.end(); ) {
    (*i)->wait();
    auto t = (*i)->compute();
    acctasks.insert(acctasks.end(), t.begin(), t.end());
    i = tasks.erase(i);

    for (auto j = acctasks.begin(); j != acctasks.end(); )
      j = (*j)->test() ? acctasks.erase(j) : ++j;
  }

  for (auto i = acctasks.begin(); i != acctasks.end(); ) {
    (*i)->wait();
    i = acctasks.erase(i); // this deallocates buffer memory
  }
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
                                shared_ptr<const Determinants> base_det, shared_ptr<const Determinants> int_det) const {

  const int norb = cc->det()->norb();
  const size_t lb = sigma->lenb();
  const size_t la = sigma->asize();

  // (astart:aend, b)
  unique_ptr<double[]> source(new double[la*lb]);
  blas::transpose(cc->local_data(), lb, la, source.get());

  unique_ptr<double[]> target(new double[la*lb]);
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

  const size_t neleb = base_det->neleb();

  const static Comb comb;
  const size_t lengb = comb(norb, neleb-2);
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
  for (auto& b : int_det->string_bits_b())
    tasks.emplace_back(la, source.get(), target.get(), hamil1.get(), base_det, b, &localmutex);

  tasks.compute();

  blas::transpose(target.get(), la, lb, source.get());
  sigma->accumulate_buffer(1.0, source);

}
