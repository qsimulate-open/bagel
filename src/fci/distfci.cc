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

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <src/fci/distfci.h>
#include <src/util/davidson.h>
#include <src/util/constants.h>
#include <vector>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace chrono;
using namespace bagel;

DistFCI::DistFCI(const multimap<string, string> a, shared_ptr<const Reference> b, const int ncore, const int nocc, const int nstate)
 : HarrisonZarrabian(a, b, ncore, nocc, nstate) {

  cout << "    * Parallel algorithm will be used." << endl;

}
    

shared_ptr<Dvec> DistFCI::form_sigma(shared_ptr<const Dvec> ccvec, shared_ptr<const MOFile> jop, const vector<int>& conv) const {
  const int ij = norb_*norb_; 

  const int nstate = ccvec->ij();

  shared_ptr<Dvec> sigmavec(new Dvec(ccvec->det(), nstate));
  sigmavec->zero();

  shared_ptr<Determinants> base_det = space_->finddet(0,0);
  shared_ptr<Determinants> int_det = space_->finddet(-1,-1);

  /* d and e are only used in the alpha-beta case and exist in the (nalpha-1)(nbeta-1) spaces */
  shared_ptr<Dvec> d(new Dvec(int_det, ij));
  shared_ptr<Dvec> e(new Dvec(int_det, ij));

  for (int istate = 0; istate != nstate; ++istate) {
    if (conv[istate]) continue;
    shared_ptr<const Civec> cc = ccvec->data(istate);  
    shared_ptr<Civec> sigma = sigmavec->data(istate);  

    vector<pair<string, double> > timing;
    auto start = high_resolution_clock::now();

    // (task1) one-electron alpha: sigma(Psib, Psi'a) += sign h'(ij) C(Psib, Psia) 
    sigma_1(cc, sigma, jop);
    if (tprint) print_timing_("task1", start, timing);

    // (task2) two electron contributions

    // (2aa) alpha-alpha contributions
    sigma_2aa(cc,sigma,jop);
    if (tprint) print_timing_("task2aa", start, timing);
    
    // (2bb) beta-beta contributions
    /* Mostly the same as the alpha-alpha, except for data storage */
    sigma_2bb(cc, sigma, jop);
    if (tprint) print_timing_("task2bb", start, timing);

    // (2ab) alpha-beta contributions
    /* Resembles more the Knowles & Handy FCI terms */
    d->zero();

    sigma_2ab(cc, sigma, jop);
    if (tprint) print_timing_("task2ab", start, timing);
    
    // (task3) one-electron beta: sigma(Psib', Psia) += sign h'(ij) C(Psib, Psia)
    sigma_3(cc, sigma, jop);

    if (tprint) {
      print_timing_("task3", start, timing);
      cout << "     timing info" << endl;
      for (auto& iter : timing )
        cout << "    " << setw(10) << iter.first << setw(10) << setprecision(2) << iter.second << endl;
    }
  }

  return sigmavec;
}


void DistFCI::sigma_2ab(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
  const int norb = norb_;

  shared_ptr<Determinants> int_det = space_->finddet(-1,-1);
  shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const int lbt = int_det->lenb();
  const int lbs = base_det->lenb();
  const double* source_base = cc->data();
  const int ij = norb*norb;

  // this is going to be a parallel loop
  for (size_t a = 0; a != base_det->lena(); ++a) {

    unique_ptr<double[]>  buf(new double[lbt*ij]);
    unique_ptr<double[]>  buf2(new double[lbt*ij]);
    fill_n(buf.get(), lbt*ij, 0.0);

    // TODO awful code
    vector<vector<DetMap>::const_iterator> aiterlist;
    for (int k = 0; k < norb; ++k)
      aiterlist.push_back(find_if(int_det->phiupa(k).begin(), int_det->phiupa(k).end(), [&a](const DetMap& i){ return i.source == a; }));

    for (int k = 0, kl = 0; k < norb; ++k) {
      // look for the element whose third element is a 
      vector<DetMap>::const_iterator aiter = aiterlist[k];

      if (aiter != int_det->phiupa(k).end()) {
        const double *source = cc->data() + aiter->target*lbs;
        for (int l = 0; l < norb; ++l, ++kl) {
          for (auto& b : int_det->phiupb(l)) {
            const double sign = aiter->sign * b.sign;
            buf[b.source+lbt*kl] += sign * source[b.target];
          }
        }
      } else {
        kl += norb;
      }
    }
    dgemm_("n", "n", lbt, ij, ij, 1.0, buf.get(), lbt, jop->mo2e_ptr(), ij, 0.0, buf2.get(), lbt);

    for (int i = 0, kl = 0; i < norb; ++i) {
      vector<DetMap>::const_iterator aiter = aiterlist[i];

      if (aiter != int_det->phiupa(i).end()) {
        double *target = sigma->data() + aiter->target*lbs;
        for (int j = 0; j < norb; ++j, ++kl) {
          for (auto& b : int_det->phiupb(j)) {
            const double sign = aiter->sign * b.sign;
            target[b.target] += sign * buf2[b.source+lbt*kl];
          }
        } 
      } else {
        kl += norb;
      }
    }
  }
}
