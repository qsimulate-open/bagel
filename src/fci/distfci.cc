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
#include <src/parallel/request.h>
#include <src/util/combination.hpp>
#include <src/util/comb.h>
#include <src/fci/hzdenomtask.h>
#include <vector>
#include <config.h>

using namespace std;
using namespace bagel;

DistFCI::DistFCI(const multimap<string, string> a, shared_ptr<const Reference> b, const int ncore, const int nocc, const int nstate)
 : FCI(a, b, ncore, nocc, nstate) {

#ifndef HAVE_MPI_H
  throw logic_error("DistFCI can be used only with MPI");
#endif

  cout << "    * Parallel algorithm will be used." << endl;

  space_ = std::shared_ptr<Space>(new Space(det_, 1));
  update(ref_->coeff());

}


vector<shared_ptr<DistCivec> > DistFCI::form_sigma(vector<shared_ptr<DistCivec> >& ccvec, shared_ptr<const MOFile> jop, const vector<int>& conv) const {
  const int ij = norb_*norb_;

  const int nstate = ccvec.size();

  vector<shared_ptr<DistCivec> > sigmavec;

  for (int istate = 0; istate != nstate; ++istate) {
    if (conv[istate]) continue;
    shared_ptr<const DistCivec> cc = ccvec[istate];
    shared_ptr<DistCivec> sigma = cc->clone(); 
    sigma->zero();

    // TODO. I don't know why we cannot initialize the queue in the constructor
    sigma->init_accum();

    vector<pair<string, double> > timing;
    Timer fcitime(1);

    sigma_aa(cc,sigma,jop);
    fcitime.tick_print("alpha-alpha");

    sigma_ab(cc, sigma, jop);
    fcitime.tick_print("alpha-beta");

    sigma_bb(cc, sigma, jop);
    fcitime.tick_print("beta-beta");

    sigma->wait();
    fcitime.tick_print("wait");

    sigmavec.push_back(sigma);
  }

  return sigmavec;
}


// accumulate Nci
void DistFCI::sigma_aa(shared_ptr<const DistCivec> cc, shared_ptr<DistCivec> sigma, shared_ptr<const MOFile> jop) const {

  shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const size_t lb = sigma->lenb();

  for (size_t a = 0; a != base_det->lena(); ++a) {
    unique_ptr<double[]> buf(new double[lb]);
    fill_n(buf.get(), lb, 0.0);
    const bitset<nbit__> nstring = base_det->stringa(a);

    for (int i = 0; i != norb_; ++i) {
      if (!nstring[i]) continue;
      bitset<nbit__> string_i = nstring; string_i.reset(i);

      // one-body term
      for (int l = 0; l != norb_; ++l) {
        if (string_i[l]) continue;
        bitset<nbit__> string_il = string_i; string_il.set(l);
        const int aloc = base_det->lexical<0>(string_il) - cc->astart();
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

            const int aloc = base_det->lexical<0>(string_ijkl) - cc->astart();
            if (aloc < 0 || aloc >= cc->asize()) continue;

            const double fac = phase * ( jop->mo2e_hz(i,j,k,l) - jop->mo2e_hz(i,j,l,k) );
            daxpy_(lb, fac, &cc->local(aloc*lb), 1, buf.get(), 1);
          }
        }
      }
    }
    // TODO this communication pattern might not be optimal
    sigma->accumulate_bstring_buf(buf, a);
    sigma->flush();
  }
}


void DistFCI::sigma_ab(shared_ptr<const DistCivec> cc, shared_ptr<DistCivec> sigma, shared_ptr<const MOFile> jop) const {

  cc->open_window();

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

    // TODO buffer should be nelec size (not norb size)
    unique_ptr<double[]>  buf2(new double[lbt*norb_*norb_]);
    unique_ptr<double[]>  buf3(new double[lbt*norb_*norb_]);
    fill_n(buf2.get(), lbt*norb_*norb_, 0.0);

    // TODO get rid of fence
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
      sigma->accumulate_bstring_buf(bcolumn, base_det->lexical<0>(atarget));
    }
    sigma->flush();
  }

  cc->close_window();
}


// beta-beta block has no communication (and should be cheap)
void DistFCI::sigma_bb(shared_ptr<const DistCivec> cc, shared_ptr<DistCivec> sigma, shared_ptr<const MOFile> jop) const {

  const shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const size_t lb = sigma->lenb();
  const size_t la = sigma->asize();

  unique_ptr<double[]> target(new double[la*lb]);
  unique_ptr<double[]> source(new double[la*lb]);

  // (astart:aend, b)
  mytranspose_(cc->local(), lb, la, source.get());
  fill_n(target.get(), la*lb, 0.0);

  // preparing Hamiltonian
  const size_t npack = norb_*(norb_-1)/2;
  unique_ptr<double[]> hamil(new double[npack*npack]);
  for (int i = 0, ijkl = 0; i != norb_; ++i)
    for (int j = 0; j < i; ++j)
      for (int k = 0; k != norb_; ++k)
        for (int l = 0; l < k; ++l, ++ijkl)
          hamil[ijkl] = jop->mo2e_hz(i,j,k,l) - jop->mo2e_hz(i,j,l,k);

  const static Comb comb;
  const size_t lengb = comb.c(norb_, neleb_-2);
  vector<bitset<nbit__> > intb(lengb, bitset<nbit__>(0));
  vector<int> data(norb_);
  iota(data.begin(), data.end(), 0);
  auto sa = intb.begin();
  do {
    for (int i=0; i < neleb_-2; ++i) sa->set(data[i]);
    ++sa;
  } while (boost::next_combination(data.begin(), data.begin()+neleb_-2, data.end()));

  // loop over intermediate string
  for (auto& b : intb) {
    // TODO buffer should be nelec size (not norb size)
    unique_ptr<double[]> ints(new double[npack*la]);
    unique_ptr<double[]> ints2(new double[npack*la]);
    fill_n(ints.get(), npack*la, 0.0);

    // first gather elements with correct sign
    for (int i = 0, ij = 0; i != norb_; ++i) {
      if (b[i]) { ij += i; continue; }
      for (int j = 0; j < i; ++j, ++ij) {
        if(b[j]) continue;
        bitset<nbit__> bt = b; bt.set(i); bt.set(j);
        const double ij_phase = base_det->sign(bt,i,j);
        daxpy_(la, ij_phase, source.get()+la*base_det->lexical<1>(bt), 1, ints.get()+la*ij, 1);
      }
    }

    // call dgemm
    dgemm_("N", "N", la, npack, npack, 1.0, ints.get(), la, hamil.get(), npack, 0.0, ints2.get(), la);

    for (int k = 0, kl = 0; k != norb_; ++k) {
      if (b[k]) { kl += k; continue; }
      for (int l = 0; l < k; ++l, ++kl) {
        if (b[l]) continue;
        bitset<nbit__> bt = b; bt.set(k); bt.set(l);
        const double kl_phase = base_det->sign(bt,k,l);
        daxpy_(la, kl_phase, ints2.get()+la*kl, 1, target.get()+la*base_det->lexical<1>(bt), 1);
      }
    }
  }
  mytranspose_(target.get(), la, lb, source.get());
  daxpy_(la*lb, 1.0, source.get(), 1, sigma->local(), 1);

  for (int i = 0; i < cc->asize(); ++i)
    for (int ip = 0; ip != nij(); ++ip)
      for (auto& iter : cc->det()->phib(ip))
        sigma->local(iter.target+cc->lenb()*i) += jop->mo1e(ip) * iter.sign * cc->local(iter.source+cc->lenb()*i);

}



void DistFCI::compute() {
  Timer pdebug(2);

  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("FCI: C1 only at the moment."); 

  // some constants
  //const int ij = nij(); 

  // Creating an initial CI vector
  vector<shared_ptr<DistCivec> > cc(nstate_);
  for (auto& i : cc)
    i = shared_ptr<DistCivec>(new DistCivec(det_));

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
    vector<shared_ptr<DistCivec> > sigma = form_sigma(cc, jop_, conv);
    pdebug.tick_print("sigma vector");

    // constructing Dvec's for Davidson
    vector<shared_ptr<const DistCivec> > ccn, sigman;
    for (auto& i : cc) ccn.push_back(i);
    for (auto& i : sigma) sigman.push_back(i);
    const vector<double> energies = davidson.compute(ccn, sigman);

    // get residual and new vectors
    vector<shared_ptr<DistCivec> > errvec = davidson.residual();
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
        if (conv[ist]) continue;
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
        davidson.orthog(c);
        list<shared_ptr<const DistCivec> > tmp;
        for (int jst = 0; jst != ist; ++jst) tmp.push_back(cc.at(jst)); 
        c->orthog(tmp);

        cc.push_back(c);
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
//cc_ = davidson.civec();
//s->print();
}


void DistFCI::update(shared_ptr<const Coeff> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  jop_ = shared_ptr<MOFile>(new Jop(ref_, ncore_, ncore_+norb_, c, "HZ"));

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
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*jop_->mo2e_hz(j, i, j, i);
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*jop_->mo2e_hz(j, i, i, j);
    }
    h[i] = jop_->mo1e(i,i);
  }
  denom_t.tick_print("jop, kop");

  denom_ = shared_ptr<DistCivec>(new DistCivec(det_));

  double* iter = denom_->local();
  vector<HZDenomTask> tasks;
  tasks.reserve(denom_->asize());
  for (size_t i = denom_->astart(); i != denom_->aend(); ++i) {
    tasks.push_back(HZDenomTask(iter, denom_->det()->stringa(i), det_, jop.get(), kop.get(), h.get()));
    iter += det()->stringb().size();
  }

  TaskQueue<HZDenomTask> tq(tasks);
  tq.compute(resources__->max_num_threads());
  denom_t.tick_print("denom");
}



vector<pair<bitset<nbit__> , bitset<nbit__> > > DistFCI::detseeds(const int ndet) {
  multimap<double, pair<size_t, size_t> > tmp;
  for (int i = 0; i != ndet; ++i)
    tmp.insert(make_pair(-1.0e10*(1+i), make_pair(0,0)));

  double* diter = denom_->local();
  for (size_t ia = denom_->astart(); ia != denom_->aend(); ++ia) {
    for (auto& b : det()->stringb()) {
      const double din = -(*diter);
      if (tmp.begin()->first < din) {
        const size_t ib = det_->lexical<1>(b);
        tmp.insert(make_pair(din, make_pair(ib, ia)));
        tmp.erase(tmp.begin());
      }
      ++diter;
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
  mpi__->allgather(&aarray[0], ndet, &aall[0], aall.size()); 
  mpi__->allgather(&barray[0], ndet, &ball[0], ball.size()); 
  mpi__->allgather(&en[0],     ndet, &eall[0], eall.size()); 

  tmp.clear();
  for (int i = 0; i != aall.size(); ++i) {
    tmp.insert(make_pair(eall[i], make_pair(ball[i], aall[i]))); 
  }

  vector<pair<bitset<nbit__> , bitset<nbit__> > > out;
  auto iter = tmp.rbegin();
  for (int i = 0; i != ndet; ++i) {
    out.push_back(make_pair(det_->stringb(iter->second.first), det_->stringa(iter->second.second)));
    ++iter;
  }

  return out;
}



// generate initial vectors
//   - bits: bit patterns of low-energy determinants
//   - nspin: #alpha - #beta
//   - out:
void DistFCI::generate_guess(const int nspin, const int nstate, vector<shared_ptr<DistCivec> > out) {
  int ndet = nstate_*10;
  start_over:
  vector<pair<bitset<nbit__>, bitset<nbit__> > > bits = detseeds(ndet);

  // Spin adapt detseeds
  int oindex = 0;
  vector<bitset<nbit__> > done;
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

    pair<vector<tuple<int, int, int> >, double> adapt = det()->spin_adapt(nelea_-neleb_, alpha, beta);
    const double fac = adapt.second;

    out[oindex]->zero();
    for (auto& ad : adapt.first) {
      const int aloc = get<1>(ad) - out[oindex]->astart();
      if (aloc >= 0 && aloc < out[oindex]->aend())
        out[oindex]->local(get<0>(ad) + det_->lenb()*aloc) = get<2>(ad)*fac;
    }

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
