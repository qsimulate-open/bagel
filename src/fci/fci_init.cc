//
// BAGEL - Parallel electron correlation program.
// Filename: fci_init.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <iomanip>
#include <stdexcept>
#include <src/fci/fci.h>
#include <src/rysint/eribatch.h>
#include <src/util/combination.hpp>
#include <iostream>

using namespace std;
using namespace bagel;

//
// generate initial vectors
//   - bits: bit patterns of low-energy determinants
//   - nspin: #alpha - #beta
//   - out:

void FCI::generate_guess(const int nspin, const int nstate, std::shared_ptr<Dvec> out) {

  int ndet = nstate_*10;
  start_over:
  vector<pair<int, int> > bits = detseeds(ndet);

  // Spin adapt detseeds
  if (true) {
    // in this case, easy. The singlet combinations are made for open-shell singlet bits
    int oindex = 0;
    vector<int> done;
    for (auto iter = bits.begin(); iter != bits.end(); ++iter) {
      const int alpha = iter->second;
      const int beta = iter->first;
      const int open_bit = (alpha^beta);

      // make sure that we have enough unpaired alpha
      const int unpairalpha = numofbits(alpha ^ (alpha & beta));
      const int unpairbeta = numofbits(beta ^ (alpha & beta));
      if (unpairalpha-unpairbeta < nelea_-neleb_) continue; 

      // check if this orbital configuration is already used
      if (find(done.begin(), done.end(), open_bit) != done.end()) continue;
      done.push_back(open_bit);

      pair<vector<tuple<int, int, int> >, double> adapt = det()->spin_adapt(nelea_-neleb_, alpha, beta);
      const double fac = adapt.second;
      for (auto iter = adapt.first.begin(); iter != adapt.first.end(); ++iter) {
        out->data(oindex)->element(get<0>(*iter), get<1>(*iter)) = get<2>(*iter)*fac;
      }

      cout << "     guess " << setw(3) << oindex << ":   closed " <<
            setw(20) << left << det()->print_bit(alpha&beta) << " open " << setw(20) << det()->print_bit(open_bit) << right << endl;

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

}


//
// returns seed determinants for initial guess
//
vector<pair<int, int> > FCI::detseeds(const int ndet) {
  multimap<double, pair<int,int> > tmp;
  for (int i = 0; i != ndet; ++i) tmp.insert(make_pair(-1.0e10*(1+i), make_pair(0,0)));

  double* diter = denom_->data();
  for (auto aiter = det()->stringa().begin(); aiter != det()->stringa().end(); ++aiter) {
    for (auto biter = det()->stringb().begin(); biter != det()->stringb().end(); ++biter, ++diter) {
      const double din = -(*diter);
      if (tmp.begin()->first < din) {
        tmp.insert(make_pair(din, make_pair(*biter, *aiter)));
        tmp.erase(tmp.begin());
      }
    }
  }
  assert(tmp.size() == ndet || ndet > det()->stringa().size()*det()->stringb().size());
  vector<pair<int, int> > out;
  for (auto iter = tmp.rbegin(); iter != tmp.rend(); ++iter)
    out.push_back(iter->second);
  return out;
}

//
// averaged diagonal elements as defined in Knowles & Handy (1989) Compt. Phys. Comm.
//
void FCI::const_denom() {

  vector<double> jop, kop, fk;
  jop.resize(norb_*norb_);
  kop.resize(norb_*norb_);
  fk.resize(norb_);
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*jop_->mo2e(j, j, i, i);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*jop_->mo2e(j, i, j, i);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    fk[i] = 0.0;
    for (int j = 0; j != norb_; ++j) {
      fk[i] += kop[i*norb_+j];
    }
  }

  shared_ptr<Civec> tmp(new Civec(det()));
  denom_ = tmp;
  const int nspin = det()->nspin();
  const int nspin2 = nspin*nspin;

  double* iter = denom_->data();
  for (auto ia = det()->stringa().begin(); ia != det()->stringa().end(); ++ia) {
    for (auto ib = det()->stringb().begin(); ib != det()->stringb().end(); ++ib, ++iter) {
      unsigned int iabit1 = *ia;
      unsigned int ibbit1 = *ib;
      const int nopen = numofbits(iabit1^ibbit1);
      const double F = (nopen >> 1) ? (static_cast<double>(nspin2 - nopen)/(nopen*(nopen-1))) : 0.0;
      *iter = 0.0;
      for (int i = 0; i != norb_; ++i, (iabit1 >>= 1), (ibbit1 >>= 1)) {
        const int nia = (iabit1&1);
        const int nib = (ibbit1&1);
        const int niab = nia + nib;
        const int Ni = (nia ^ nib);
        unsigned int iabit2 = *ia;
        unsigned int ibbit2 = *ib;
        for (int j = 0; j != i; ++j, (iabit2 >>= 1), (ibbit2 >>= 1)) {
          const int nja = (iabit2&1);
          const int njb = (ibbit2&1);
          const int Nj = (nja ^ njb);
          const int addj = niab * (nja + njb);
          *iter += jop[j+norb_*i] * 2.0 * addj - kop[j+norb_*i] * (F*Ni*Nj + addj);
        }
        *iter += (jop_->mo1e(i,i) + fk[i]) * niab - kop[i+norb_*i] * 0.5 * (Ni - niab*niab);
      }
    }
  }
}


void FCI::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "        FCI calculation      " << endl;
  cout << "  ---------------------------" << endl << endl;
}

