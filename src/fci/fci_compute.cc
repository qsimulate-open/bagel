//
// Newint - Parallel electron correlation program.
// Filename: fci_compute.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <src/fci/fci.h>
#include <src/util/davidson.h>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;

static const int unit = 1;
static const double one = 1.0;
static const double half = 0.5;
static const double zero = 0.0; 
static const string indent = "  ";
static const string space3 = "   "; 

static int address(int i, int j) { assert(i <= j); return i+((j*(j+1))>>1); };

void FCI::compute() {

  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("FCI: C1 only at the moment."); 

  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  jop_ = shared_ptr<MOFile>(new MOFile(geom_, ref_));

  // right now full basis is used. 
  int start_int = ::clock();
  core_energy_ = jop_->create_Jiiii(ncore_, ncore_+norb_, ncore_);
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) <<
          static_cast<double>(::clock() - start_int)/static_cast<double>(CLOCKS_PER_SEC) << endl << endl;

  // create denominator here. Stored in shared<Civec> denom_
  const_denom();

  // some constants
  const int la = stringa_.size();
  const int lb = stringb_.size();
  const int ij = norb_ * (norb_+1) /2;

  // we need two vectors for intermediate quantities
  shared_ptr<Dvec> d(new Dvec(lb, la, ij));
  shared_ptr<Dvec> e(new Dvec(lb, la, ij));

  // Creating an initial CI vector
  shared_ptr<Dvec> cc_tmp(new Dvec(lb, la, nstate_)); // B runs first
  cc_ = cc_tmp;
  shared_ptr<Dvec> sigma(new Dvec(lb, la, nstate_));

  // find determinants that have small diagonal energies
  generate_guess(nelea_-neleb_, nstate_, cc_); 
  // TODO note that generate_guess is only working fine for singlets

  // nuclear energy retrieved from geometry
  const double nuc_core = geom_->nuclear_repulsion() + core_energy();

  // Davidson utility
  DavidsonDiag<Civec> davidson(nstate_, max_iter_);

  // main iteration starts here
  cout << "  === FCI iteration ===" << endl << endl;
  // 0 means not converged
  vector<int> conv(nstate_,0);

  for (int iter = 0; iter != max_iter_; ++iter) { 
    int start = ::clock();

    // form a sigma vector given cc
    form_sigma(cc_, sigma, d, e, conv);

    // constructing Dvec's for Davidson
    shared_ptr<const Dvec> ccn(new Dvec(cc_));
    shared_ptr<const Dvec> sigman(new Dvec(sigma));
    const vector<double> energies = davidson.compute(ccn->dvec(conv), sigman->dvec(conv));

    // get residual and new vectors
    vector<shared_ptr<Civec> > errvec = davidson.residual();

    // compute errors
    vector<double> errors;
    for (int i = 0; i != nstate_; ++i) {
      errors.push_back(errvec[i]->variance());
      conv[i] = static_cast<int>(errors[i] < thresh_);
    }

    if (!*min_element(conv.begin(), conv.end())) {
      // denominator scaling 
      for (int ist = 0; ist != nstate_; ++ist) {
        if (conv[ist]) continue;
        const int size = cc_->data(ist)->size();
        double* target_array = cc_->data(ist)->first();
        double* source_array = errvec[ist]->first();
        double* denom_array = denom_->first();
        const double en = energies[ist];
        for (int i = 0; i != size; ++i) {
          target_array[i] = source_array[i] / min(en - denom_array[i], -0.1);
        }
        davidson.orthog(cc_->data(ist));
        list<shared_ptr<const Civec> > tmp;
        for (int jst = 0; jst != ist; ++jst) tmp.push_back(cc_->data(jst)); 
        cc_->data(ist)->orthog(tmp);
      }
    }

    // printing out
    int end = ::clock();
    if (nstate_ != 1 && iter) cout << endl;
    for (int i = 0; i != nstate_; ++i) {
      cout << indent << setw(5) << iter << setw(3) << i << setw(2) << (conv[i] ? "*" : " ")
                                        << setw(17) << fixed << setprecision(8) << energies[i]+nuc_core << space3 
                                        << setw(10) << scientific << setprecision(2) << errors[i] << fixed << setw(10) << setprecision(2)
                                        << (end - start)/static_cast<double>(CLOCKS_PER_SEC) << endl; 
      energy_[i] = energies[i]+nuc_core;
    }
    if (*min_element(conv.begin(), conv.end())) break;
  }
  // main iteration ends here

  vector<shared_ptr<Civec> > s = davidson.civec();
  print_civectors(s, 0.05);
  cc_ = shared_ptr<Dvec>(new Dvec(s));

}




void FCI::form_sigma(shared_ptr<Dvec> ccvec, shared_ptr<Dvec> sigmavec,
                     shared_ptr<Dvec> d, shared_ptr<Dvec> e,
                     const vector<int>& conv) { // d and e are scratch area for D and E intermediates 
  const int la = d->lena();  
  const int lb = d->lenb();  
  const int ij = d->ij();
  sigmavec->zero();

  const int nstate = ccvec->ij();

  for (int istate = 0; istate != nstate; ++istate) {
    if (conv[istate]) continue;
    shared_ptr<Civec> cc = ccvec->data(istate);  
    shared_ptr<Civec> sigma = sigmavec->data(istate);  

    vector<pair<string, double> > timing;
    int start = ::clock();

    // (task1) one-electron alpha: sigma(Psib, Psi'a) += sign h'(ij) C(Psib, Psia) 
    sigma_1(cc, sigma);
    if (tprint) print_timing_("task1", start, timing);

    // (task2) two electron contributions
    d->zero();

    // step (c) (task2a-1) D(Phib, Phia, ij) += sign C(Psib, Phi'a)
    sigma_2a1(cc, d);
    if (tprint) print_timing_("task2a-1", start, timing);

    // step (d) (task2a-2) D(Phib, Phia, ij) += sign C(Psib', Phia)
    sigma_2a2(cc, d);
    if (tprint) print_timing_("task2a-2", start, timing);

    // step (e) (task2b) E(Phib, Phia, ij) = D(Psib, Phia, ij) (ij|kl)
    sigma_2b(d, e);
    if (tprint) print_timing_("task2b", start, timing);

    // step (f) (task2c-1) sigma(Phib, Phia') += sign E(Psib, Phia, ij)
    sigma_2c1(sigma, e);
    if (tprint) print_timing_("task2c-1", start, timing);

    // step (g) (task2c-2) sigma(Phib', Phia) += sign E(Psib, Phia, ij)
    sigma_2c2(sigma, e);
    if (tprint) print_timing_("task2c-2", start, timing);

    // (task3) one-electron beta: sigma(Psib', Psia) += sign h'(ij) C(Psib, Psia)
    sigma_3(cc, sigma);

    if (tprint) {
      print_timing_("task3", start, timing);
      cout << "     timing info" << endl;
      for (auto iter = timing.begin(); iter != timing.end(); ++iter)
        cout << "    " << setw(10) << iter->first << setw(10) << setprecision(2) << iter ->second << endl;
    }
  }
}


void FCI::sigma_1(shared_ptr<Civec> cc, shared_ptr<Civec> sigma) {
  const int ij = jop_->sizeij();
  const int lb = cc->lenb();
  for (int ip = 0; ip != ij; ++ip) {
    const double h = jop_->mo1e(ip);
    for (auto iter = phia_[ip].begin();  iter != phia_[ip].end(); ++iter) {
      const double hc = h * get<1>(*iter);
      daxpy_(&lb, &hc, cc->element_ptr(0, get<2>(*iter)), &unit, sigma->element_ptr(0, get<0>(*iter)), &unit); 
    }
  }
}

void FCI::sigma_2a1(shared_ptr<Civec> cc, shared_ptr<Dvec> d) {
  const int lb = d->lenb();
  const int ij = d->ij();
  double* const source_base = cc->first();
  for (int ip = 0; ip != ij; ++ip) {
    double* const target_base = d->data(ip)->first();
    for (auto iter = phia_[ip].begin();  iter != phia_[ip].end(); ++iter) {
      const double sign = static_cast<double>(get<1>(*iter));
      double* const target_array = target_base + get<2>(*iter)*lb;
      daxpy_(&lb, &sign, source_base + get<0>(*iter)*lb, &unit, target_array, &unit);
    }
  }
}

void FCI::sigma_2a2(shared_ptr<Civec> cc, shared_ptr<Dvec> d) {
  const int la = d->lena();
  const int ij = d->ij();
  for (int i = 0; i < la; i+=4) {
    if (i+3 < la) {
      double* const source_array0 = cc->element_ptr(0, i); 
      double* const source_array1 = cc->element_ptr(0, i+1); 
      double* const source_array2 = cc->element_ptr(0, i+2); 
      double* const source_array3 = cc->element_ptr(0, i+3); 
      for (int ip = 0; ip != ij; ++ip) {
        double* const target_array0 = d->data(ip)->element_ptr(0, i); 
        double* const target_array1 = d->data(ip)->element_ptr(0, i+1); 
        double* const target_array2 = d->data(ip)->element_ptr(0, i+2); 
        double* const target_array3 = d->data(ip)->element_ptr(0, i+3); 
        for (auto iter = phib_[ip].begin(); iter != phib_[ip].end(); ++iter) {
          const double sign = static_cast<double>(get<1>(*iter));
          const int ia = get<2>(*iter);
          const int ib = get<0>(*iter);
          target_array0[ia] += sign * source_array0[ib]; 
          target_array1[ia] += sign * source_array1[ib]; 
          target_array2[ia] += sign * source_array2[ib]; 
          target_array3[ia] += sign * source_array3[ib]; 
        }
      } 
    } else {
      for (int j = i; j != la; ++j) {
        double* const source_array0 = cc->element_ptr(0, j);
        for (int ip = 0; ip != ij; ++ip) {
          double* const target_array0 = d->data(ip)->element_ptr(0, j);
          for (auto iter = phib_[ip].begin(); iter != phib_[ip].end(); ++iter) {
            const double sign = static_cast<double>(get<1>(*iter));
            target_array0[get<2>(*iter)] += sign * source_array0[get<0>(*iter)]; 
          }
        }
      }
    }
  }
}

void FCI::sigma_2c1(shared_ptr<Civec> sigma, shared_ptr<Dvec> e) {
  const int lb = e->lenb();
  const int ij = e->ij();
  for (int ip = 0; ip != ij; ++ip) { 
    double* const source_base = e->data(ip)->first();
    for (auto iter = phia_[ip].begin(); iter != phia_[ip].end(); ++iter) {
      const double sign = static_cast<double>(get<1>(*iter)); 
      double* const target_array = sigma->element_ptr(0, get<0>(*iter));
      daxpy_(&lb, &sign, source_base + lb*get<2>(*iter), &unit, target_array, &unit);
    }
  }
}

void FCI::sigma_2c2(shared_ptr<Civec> sigma, shared_ptr<Dvec> e) {
  const int la = e->lena();
  const int ij = e->ij();
  for (int i = 0; i < la; i+=8) {
    if (i+7 < la) {
      double* const target_array0 = sigma->element_ptr(0, i);
      double* const target_array1 = sigma->element_ptr(0, i+1);
      double* const target_array2 = sigma->element_ptr(0, i+2);
      double* const target_array3 = sigma->element_ptr(0, i+3);
      double* const target_array4 = sigma->element_ptr(0, i+4);
      double* const target_array5 = sigma->element_ptr(0, i+5);
      double* const target_array6 = sigma->element_ptr(0, i+6);
      double* const target_array7 = sigma->element_ptr(0, i+7);
      for (int ip = 0; ip != ij; ++ip) {
        double* const source_array0 = e->data(ip)->element_ptr(0, i);
        double* const source_array1 = e->data(ip)->element_ptr(0, i+1);
        double* const source_array2 = e->data(ip)->element_ptr(0, i+2);
        double* const source_array3 = e->data(ip)->element_ptr(0, i+3);
        double* const source_array4 = e->data(ip)->element_ptr(0, i+4);
        double* const source_array5 = e->data(ip)->element_ptr(0, i+5);
        double* const source_array6 = e->data(ip)->element_ptr(0, i+6);
        double* const source_array7 = e->data(ip)->element_ptr(0, i+7);
        for (auto iter = phib_[ip].begin(); iter != phib_[ip].end(); ++iter) {
          const double sign = static_cast<double>(get<1>(*iter));
          target_array0[get<0>(*iter)] += sign * source_array0[get<2>(*iter)]; 
          target_array1[get<0>(*iter)] += sign * source_array1[get<2>(*iter)]; 
          target_array2[get<0>(*iter)] += sign * source_array2[get<2>(*iter)]; 
          target_array3[get<0>(*iter)] += sign * source_array3[get<2>(*iter)]; 
          target_array4[get<0>(*iter)] += sign * source_array4[get<2>(*iter)]; 
          target_array5[get<0>(*iter)] += sign * source_array5[get<2>(*iter)]; 
          target_array6[get<0>(*iter)] += sign * source_array6[get<2>(*iter)]; 
          target_array7[get<0>(*iter)] += sign * source_array7[get<2>(*iter)]; 
        }
      }
    } else {
      for (int j = i; j != la; ++j) {
        double* const target_array0 = sigma->element_ptr(0, j);
        for (int ip = 0; ip != ij; ++ip) {
          double* const source_array0 = e->data(ip)->element_ptr(0, j);
          for (auto iter = phib_[ip].begin(); iter != phib_[ip].end(); ++iter) {
            const double sign = static_cast<double>(get<1>(*iter));
            target_array0[get<0>(*iter)] += sign * source_array0[get<2>(*iter)]; 
          }
        }
      }
    }
  }
}


void FCI::sigma_3(shared_ptr<Civec> cc, shared_ptr<Civec> sigma) {
  const int la = cc->lena();
  const int ij = jop_->sizeij();

  for (int i = 0; i < la; i+=8) {
    if (i+7 < la) {
      double* const target_array0 = sigma->element_ptr(0, i);
      double* const target_array1 = sigma->element_ptr(0, i+1);
      double* const target_array2 = sigma->element_ptr(0, i+2);
      double* const target_array3 = sigma->element_ptr(0, i+3);
      double* const target_array4 = sigma->element_ptr(0, i+4);
      double* const target_array5 = sigma->element_ptr(0, i+5);
      double* const target_array6 = sigma->element_ptr(0, i+6);
      double* const target_array7 = sigma->element_ptr(0, i+7);
      double* const source_array0 = cc->element_ptr(0, i);
      double* const source_array1 = cc->element_ptr(0, i+1);
      double* const source_array2 = cc->element_ptr(0, i+2);
      double* const source_array3 = cc->element_ptr(0, i+3);
      double* const source_array4 = cc->element_ptr(0, i+4);
      double* const source_array5 = cc->element_ptr(0, i+5);
      double* const source_array6 = cc->element_ptr(0, i+6);
      double* const source_array7 = cc->element_ptr(0, i+7);
      for (int ip = 0; ip != ij; ++ip) {
        const double h = jop_->mo1e(ip);
        for (auto iter = phib_[ip].begin();  iter != phib_[ip].end(); ++iter) {
          const double hc = h * get<1>(*iter);
          target_array0[get<0>(*iter)] += hc * source_array0[get<2>(*iter)];
          target_array1[get<0>(*iter)] += hc * source_array1[get<2>(*iter)];
          target_array2[get<0>(*iter)] += hc * source_array2[get<2>(*iter)];
          target_array3[get<0>(*iter)] += hc * source_array3[get<2>(*iter)];
          target_array4[get<0>(*iter)] += hc * source_array4[get<2>(*iter)];
          target_array5[get<0>(*iter)] += hc * source_array5[get<2>(*iter)];
          target_array6[get<0>(*iter)] += hc * source_array6[get<2>(*iter)];
          target_array7[get<0>(*iter)] += hc * source_array7[get<2>(*iter)];
        }
      }
    } else {
      for (int j = i; j != la; ++j) {  
        double* const target_array0 = sigma->element_ptr(0, j);
        double* const source_array0 = cc->element_ptr(0, j);
        for (int ip = 0; ip != ij; ++ip) {
          const double h = jop_->mo1e(ip);
          for (auto iter = phib_[ip].begin();  iter != phib_[ip].end(); ++iter) {
            const double hc = h * get<1>(*iter);
            target_array0[get<0>(*iter)] += hc * source_array0[get<2>(*iter)];
          }
        }
      }
    }
  }
}


void FCI::sigma_2b(shared_ptr<Dvec> d, shared_ptr<Dvec> e) {
  const int la = d->lena();
  const int lb = d->lenb();
  const int ij = d->ij();
  const int lenab = la*lb;
  dgemm_("n", "n", &lenab, &ij, &ij, &half, d->first(), &lenab, jop_->mo2e_ptr(), &ij,
                                     &zero, e->first(), &lenab);
}

