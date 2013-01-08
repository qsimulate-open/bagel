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
//#include <src/parallel/request.h>
#include <vector>
#include <config.h>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace chrono;
using namespace bagel;

DistFCI::DistFCI(const multimap<string, string> a, shared_ptr<const Reference> b, const int ncore, const int nocc, const int nstate)
 : HarrisonZarrabian(a, b, ncore, nocc, nstate) {

  cout << "    * Parallel algorithm will be used." << endl;
#ifndef HAVE_MPI_H
  throw logic_error("DistFCI can be used only with MPI");
#endif

}
    

shared_ptr<Dvec> DistFCI::form_sigma(shared_ptr<const Dvec> ccvec, shared_ptr<const MOFile> jop, const vector<int>& conv) const {
  const int ij = norb_*norb_; 

  const int nstate = ccvec->ij();

  shared_ptr<Dvec> sigmavec(new Dvec(ccvec->det(), nstate));
  sigmavec->zero();

  shared_ptr<Determinants> base_det = space_->finddet(0,0);
  shared_ptr<Determinants> int_det = space_->finddet(-1,-1);

  for (int istate = 0; istate != nstate; ++istate) {
    if (conv[istate]) continue;
    shared_ptr<const Civec> cc = ccvec->data(istate);  
    shared_ptr<Civec> sigma = sigmavec->data(istate);  

    vector<pair<string, double> > timing;
    Timer fcitime(1);

    // (1aa, 2aa) alpha-alpha contributions
    sigma_aa(cc,sigma,jop);
    fcitime.tick_print("alpha-alpha");
    
    // (1bb, 2bb) beta-beta contributions
    sigma_bb(cc, sigma, jop);
    fcitime.tick_print("beta-beta");

    // (2ab) alpha-beta contributions
    sigma_2ab(cc, sigma, jop);
    fcitime.tick_print("alpha-beta");
  }

  return sigmavec;
}


// accumulate Nci
void DistFCI::sigma_aa(shared_ptr<const Civec> ccg, shared_ptr<Civec> sigmag, shared_ptr<const MOFile> jop) const {

  shared_ptr<const DistCivec> cc = ccg->distcivec();
  shared_ptr<DistCivec> sigma = cc->clone();

  sigma->open_window();

  shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const int lb = sigma->lenb();

  unique_ptr<double[]> buf(new double[lb]);
  
  for (size_t a = 0; a != base_det->lena(); ++a) {
    fill_n(buf.get(), lb, 0.0);
    const bitset<nbit__> nstring = base_det->stringa(a);

    for (int i = 0; i != norb_; ++i) {
      if (!nstring[i]) continue;
      bitset<nbit__> string_i = nstring; string_i.reset(i);

      // one-body term
      for (int l = 0; l != norb_; ++l) {
        if (string_i[l]) continue;
        bitset<nbit__> string_il = string_i; string_il.set(l);
        const size_t aloc = base_det->lexical<0>(string_il) - cc->astart();
        if (aloc < 0 || aloc >= cc->asize()) continue; 

        const double fac = jop->mo1e(min(i,l)+max(i,l)*(max(i,l)+1)/2) * base_det->sign(string_il,l,i);
        daxpy_(lb, fac, &cc->local(aloc*lb), 1, buf.get(), 1);
      }

      // two-body term
      for (int j = 0; j < i; ++j) {
        if(!nstring[j]) continue;
        const int ij_phase = base_det->sign(nstring,i,j);
        bitset<nbit__> string_ij = nstring;
        string_ij.reset(i);
        string_ij.reset(j);
        for (int l = 0; l != norb_; ++l) {
          if (string_ij[l]) continue;
          for (int k = 0; k < l; ++k) {
            if (string_ij[k]) continue;
            const int kl_phase = base_det->sign(string_ij,l,k);
            const double phase = - (ij_phase*kl_phase);
            bitset<nbit__> string_ijkl = string_ij;
            string_ijkl.set(k); string_ijkl.set(l);

            const size_t aloc = base_det->lexical<0>(string_ijkl) - cc->astart();
            if (aloc < 0 || aloc >= cc->asize()) continue; 

            const double fac = phase * ( jop->mo2e_hz(i,j,k,l) - jop->mo2e_hz(i,j,l,k) );
            daxpy_(lb, fac, &cc->local(aloc*lb), 1, buf.get(), 1);
          }
        }
      }
    }
    sigma->accumulate_bstring(buf.get(), a);
  }
  sigma->close_window();
  *sigmag += *sigma->civec();
}


void DistFCI::sigma_2ab(shared_ptr<const Civec> ccg, shared_ptr<Civec> sigmag, shared_ptr<const MOFile> jop) const {
  shared_ptr<const DistCivec> cc = ccg->distcivec();
  shared_ptr<DistCivec> sigma = cc->clone();

  cc->open_window();
  sigma->open_window();

  shared_ptr<Determinants> int_det = space_->finddet(-1,-1);
  shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const size_t lbt = int_det->lenb();
  const size_t lbs = base_det->lenb();
  const int ij = norb_*norb_;

  const int rank = mpi__->rank();
  const int size = mpi__->size();

  const size_t nloop = (int_det->lena()-1)/size+1;

  // shamelessly statically distributing across processes
  for (size_t loop = 0; loop != nloop; ++loop) {

    size_t a = rank + loop*size;
    if (a >= int_det->lena()) { cc->fence(); break; } // fence needed, otherwise stall 

    const bitset<nbit__> astring = int_det->stringa(a);

    // first receive all the data (nele_a * lenb)
    unique_ptr<double[]>  buf(new double[lbs*norb_]);
    fill_n(buf.get(), lbs*norb_, 0.0);

    for (int i = 0; i != norb_; ++i) {
      if (!astring[i]) {
        bitset<nbit__> tmp = astring; tmp.set(i);
        cc->get_bstring(buf.get()+i*lbs, base_det->lexical<0>(tmp));
      }
    }

    // TODO can be smaller
    unique_ptr<double[]>  buf2(new double[lbt*norb_*norb_]);
    unique_ptr<double[]>  buf3(new double[lbt*norb_*norb_]);
    fill_n(buf2.get(), lbt*norb_*norb_, 0.0);

    // somehow I need to do this here, which means that all the processes sync here
    // (which might not be too bad for static load balancing)
    cc->fence();

    for (int k = 0, kl = 0; k != norb_; ++k)
      for (int l = 0; l != norb_; ++l, ++kl)
        for (auto& b : int_det->phiupb(l))
          buf2[b.source+lbt*kl] += base_det->sign(astring, -1, k) * b.sign * buf[b.target+k*lbs];

    dgemm_("n", "n", lbt, ij, ij, 1.0, buf2.get(), lbt, jop->mo2e_ptr(), ij, 0.0, buf3.get(), lbt);

    for (int i = 0; i < norb_; ++i) {
      if (astring[i]) continue;
      bitset<nbit__> atarget = astring; atarget.set(i);
      const double asign = base_det->sign(astring, -1, i);
      
      unique_ptr<double[]> bcolumn(new double[lbs]);
      fill_n(bcolumn.get(), lbs, 0.0);

      for (int j = 0; j < norb_; ++j) {
        for (auto& b : int_det->phiupb(j))
          bcolumn[b.target] += asign * b.sign * buf3[b.source+lbt*(j+norb_*i)];
      }
      sigma->accumulate_bstring(bcolumn.get(), base_det->lexical<0>(atarget));
    }
  }

  sigma->close_window();
  cc->close_window();
  *sigmag += *sigma->civec();
}


// beta-beta block has no communication (and should be cheap)
void DistFCI::sigma_bb(shared_ptr<const Civec> ccg, shared_ptr<Civec> sigmag, shared_ptr<const MOFile> jop) const {

  shared_ptr<const DistCivec> cc = ccg->distcivec();
  shared_ptr<DistCivec> sigma = cc->clone();

  const shared_ptr<Determinants> base_det = space_->finddet(0,0);
  const double* const source_base = cc->local();
  double* target_base = sigma->local();

  const int lb = sigma->lenb();
  for (auto biter = base_det->stringb().begin(); biter != base_det->stringb().end(); ++biter,++target_base) {
    bitset<nbit__> nstring = *biter;
    for (int i = 0; i != norb_; ++i) {
      if (!nstring[i]) continue;
      for (int j = 0; j < i; ++j) {
        if(!nstring[j]) continue;
        const int ij_phase = base_det->sign(nstring,i,j);
        bitset<nbit__> string_ij = nstring;
        string_ij.reset(i); string_ij.reset(j);
        for (int l = 0; l != norb_; ++l) {
          if (string_ij[l]) continue;
          for (int k = 0; k < l; ++k) {
            if (string_ij[k]) continue;
            const int kl_phase = base_det->sign(string_ij,l,k);
            const double phase = -static_cast<double>(ij_phase*kl_phase);
            bitset<nbit__> string_ijkl = string_ij;
            string_ijkl.set(k); string_ijkl.set(l);
            const double temp = phase * ( jop->mo2e_hz(i,j,k,l) - jop->mo2e_hz(i,j,l,k) );
            const double* source = source_base + base_det->lexical<1>(string_ijkl);
            daxpy_(sigma->asize(), temp, source, lb, target_base, lb);
          }
        }
      }
    }
  }
  for (int i = 0; i < cc->asize(); ++i)
    for (int ip = 0; ip != nij(); ++ip)
      for (auto& iter : cc->det()->phib(ip))
        sigma->local(iter.target+cc->lenb()*i) += jop->mo1e(ip) * iter.sign * cc->local(iter.source+cc->lenb()*i);

  *sigmag += *sigma->civec();
}

