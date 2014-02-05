//
// BAGEL - Parallel electron correlation program.
// Filename: ras/rasci_denom.cc
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
#include <src/ras/distrasci.h>

using namespace std;
using namespace bagel;

namespace bagel {
  namespace RAS {
    struct DenomTask {
      double* const data_;
      const bitset<nbit__> abit_;
      shared_ptr<const RASString> stringb_;
      const double* const jop_;
      const double* const kop_;
      const double* const h_;

      DenomTask(double* o, bitset<nbit__> ia, shared_ptr<const RASString> sb, double* j, double* k, double* h) :
        data_(o), abit_(ia), stringb_(sb), jop_(j), kop_(k), h_(h) {}

      void compute() {
        const int nspin = abit_.count() - stringb_->nele();
        const int nspin2 = nspin*nspin;
        const int norb = stringb_->norb();

        double* iter = data_;
        const bitset<nbit__> ia = abit_;
        for (auto& ib : *stringb_) {
          const int nopen = (ia^ib).count();
          const double F = (nopen >> 1) ? (static_cast<double>(nspin2 - nopen)/(nopen*(nopen-1))) : 0.0;
          *iter = 0.0;
          for (int i = 0; i < norb; ++i) {
            const int nia = ia[i];
            const int nib = ib[i];
            const int niab = (nia + nib);
            const int Ni = (nia ^ nib);
            for (int j = 0; j < i; ++j) {
              const int nja = ia[j];
              const int njb = ib[j];
              const int Nj = nja ^ njb;
              const int addj = niab * (nja + njb);

              *iter += jop_[j+norb*i] * 2.0 * addj - kop_[j+norb*i] * (F*Ni*Nj + addj);
            }
            *iter += h_[i] * niab - kop_[i+norb*i] * 0.5 * (Ni - niab*niab);
          }
          ++iter;
        }
      }
    };
  }
}

void RASCI::const_denom() {
  Timer denom_t;
  unique_ptr<double[]> h(new double[norb_]);
  unique_ptr<double[]> jop(new double[norb_*norb_]);
  unique_ptr<double[]> kop(new double[norb_*norb_]);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*jop_->mo2e_hz(j, i, j, i);
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*jop_->mo2e_hz(j, i, i, j);
    }
    h[i] = jop_->mo1e(i,i);
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<RASCivec>(det());

  size_t tasksize = 0;
  for (auto& iblock : denom_->blocks()) { if (iblock) tasksize += iblock->lena(); }
  TaskQueue<RAS::DenomTask> tasks(tasksize);

  for (auto& iblock : denom_->blocks()) {
    if ( !iblock ) continue;
    double* iter = iblock->data();
    for (auto& ia : *iblock->stringa()) {
      tasks.emplace_back(iter, ia, iblock->stringb(), jop.get(), kop.get(), h.get());
      iter += iblock->lenb();
    }
  }

  tasks.compute();
  denom_t.tick_print("denom");
}

void DistRASCI::const_denom() {
  Timer denom_t;
  unique_ptr<double[]> h(new double[norb_]);
  unique_ptr<double[]> jop(new double[norb_*norb_]);
  unique_ptr<double[]> kop(new double[norb_*norb_]);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*jop_->mo2e_hz(j, i, j, i);
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*jop_->mo2e_hz(j, i, i, j);
    }
    h[i] = jop_->mo1e(i,i);
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<DistRASCivec>(det_);

  size_t tasksize = 0;
  for (auto& iblock : denom_->blocks()) { if (iblock) tasksize += iblock->asize(); }
  TaskQueue<RAS::DenomTask> tasks(tasksize);

  for (auto& iblock : denom_->blocks()) {
    if ( !iblock ) continue;
    double* iter = iblock->local();
    for (size_t ia = iblock->astart(); ia < iblock->aend(); ++ia, iter+=iblock->lenb())
      tasks.emplace_back(iter, det_->stringa(ia + iblock->stringa()->offset()), iblock->stringb(), jop.get(), kop.get(), h.get());
  }

  tasks.compute();
  denom_t.tick_print("denom");
}

void RASCI::update(shared_ptr<const Coeff> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  // Same Jop as used in FCI
  jop_ = make_shared<Jop>(ref_, ncore_, ncore_+norb_, c, "HZ");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;

  const_denom();
}

void DistRASCI::update(shared_ptr<const Coeff> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  // Same Jop as used in FCI
  jop_ = make_shared<Jop>(ref_, ncore_, ncore_+norb_, c, "HZ");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;

  const_denom();
}
