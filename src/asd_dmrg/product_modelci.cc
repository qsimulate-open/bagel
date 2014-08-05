//
// BAGEL - Parallel electron correlation program.
// Filename: product_modelci.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <iostream>
#include <algorithm>

#include <src/asd_dmrg/product_modelci.h>
#include <src/asd_dmrg/product_citask.h>

using namespace std;
using namespace bagel;

double ProductCIHamTask::compute_pure_ras(const bitset<nbit__> abra, const bitset<nbit__> bbra, const bitset<nbit__> aket, const bitset<nbit__> bket) {
  assert(abra.count()==aket.count() && bbra.count()==bket.count());

  const bitset<nbit__> aexch = abra ^ aket;
  const bitset<nbit__> bexch = bbra ^ bket;

  const int naexch = aexch.count();
  const int nbexch = bexch.count();
  const int nexch = naexch + nbexch;

  double out = 0.0;

  if (nexch == 0) {
    // diagonal contribution
    vector<int> aoccs = bagel::bit_to_numbers(abra);
    vector<int> boccs = bagel::bit_to_numbers(bbra);

    // one-body
    for_each(aoccs.begin(), aoccs.end(), [this, &out] (const int& j) { out += mo1e_->element(j,j); });
    for_each(boccs.begin(), boccs.end(), [this, &out] (const int& j) { out += mo1e_->element(j,j); });

    // two-body
    for (auto i = aoccs.begin(); i != aoccs.end(); ++i) {
      for (auto j = aoccs.begin(); j != i; ++j)
        out += jop_->mo2e(*j,*j,*i,*i) - jop_->mo2e(*j,*i,*j,*i);
      for (auto j = boccs.begin(); j != boccs.end(); ++j)
        out += jop_->mo2e(*j,*j,*i,*i);
    }

    for (auto i = boccs.begin(); i != boccs.end(); ++i)
      for (auto j = boccs.begin(); j != i; ++j)
        out += jop_->mo2e(*j,*j,*i,*i) - jop_->mo2e(*j,*i,*j,*i);
  }
  else if (nexch == 2) {
    // single exchange
    vector<int> exch_ind = bagel::bit_to_numbers(naexch==2 ? aexch : bexch);
    const int i = exch_ind.front();
    const int j = exch_ind.back();
    const double signij = static_cast<double>( bagel::sign((naexch==2?abra:bbra), i, j) );

    // one-body
    out += signij * mo1e_->element(i,j);

    // two-body
    bitset<nbit__> acomm = abra & aket;
    bitset<nbit__> bcomm = bbra & bket;
    for (int k = 0; k < rnorb_; ++k) {
      const double nj = static_cast<double>(acomm[k] + bcomm[k]);
      out += nj * signij * jop_->mo2e(i, j, k, k);
    }

    for (int& k : bit_to_numbers((naexch==2?acomm:bcomm)))
      out -= signij * mo2e(i, k, j, k);
  }
  else if (nexch == 4) {
    // double exchange (two-body only)
    if (naexch == nbexch) {
      vector<int> alphas = bagel::bit_to_numbers(aexch);
      vector<int> betas  = bagel::bit_to_numbers(bexch);
      const int i = alphas.front();
      const int j = alphas.back();
      const int k = betas.front();
      const int l = betas.back();

      const int signij = bagel::sign(abra, i, j);
      const int signkl = bagel::sign(bbra, k, l);
      out += signij * signkl * jop_->mo2e(i,j,k,l);
    }
    else {
      bitset<nbit__> exch = (naexch==4 ? aexch : bexch);
      bitset<nbit__> eket = (naexch==4 ? aket : bket);
      bitset<nbit__> ebra = (naexch==4 ? abra : bbra);
      vector<int> ann_list = bagel::bit_to_numbers(exch & eket);
      vector<int> cre_list = bagel::bit_to_numbers(exch & ebra);
      const int i = ann_list.front();
      const int j = ann_list.back();
      const int k = cre_list.front();
      const int l = cre_list.back();
      bitset<nbit__> tmp = eket;
      tmp.reset(i); tmp.reset(j);
      const double phase = static_cast<double>(bagel::sign(eket, i, j) * bagel::sign(tmp, k, l));
      out += phase * (mo2e(i,k,j,l) - mo2e(i,l,k,j));
    }
  }
  return out;
}

double ProductCIHamTask::matrix_element_impl(const PCI::Basis& bra, const PCI::Basis& ket) {
  assert((bra.nelea+bra.neleb+bra.alpha.count()+bra.beta.count()) == (ket.nelea+ket.neleb+ket.alpha.count()+ket.beta.count()));
  if ( !(bra.alpha.count()==ket.alpha.count() && bra.beta.count()==ket.beta.count()) )
    throw logic_error("non-particle-number conserving terms not yet implemented!");

  const int ntransa = bra.alpha.count() - ket.alpha.count();
  const int ntransb = bra.beta.count() - ket.beta.count();

  double out = 0.0;

  if (ntransa==0 && ntransb==0) {// "Diagonal contribution"
    if (bra.state==ket.state)
      out += compute_pure_ras(bra.alpha, bra.beta, ket.alpha, ket.beta);

    if (make_pair(bra.alpha,bra.beta)==make_pair(ket.alpha,ket.beta)) // |L>==|L'> --> <phi| H |phi'>
      out += (*blockops_->ham(bra.key()))(bra.state, ket.state);

    shared_ptr<const btas::Tensor4<double>> Qaa = blockops_->Q_aa(bra.key());
    shared_ptr<const btas::Tensor4<double>> Qbb = blockops_->Q_bb(bra.key());

    const bitset<nbit__> aexch = bra.alpha ^ ket.alpha;
    const bitset<nbit__> bexch = bra.beta ^ ket.beta;

    const int naexch = aexch.count();
    const int nbexch = bexch.count();
    const int nexch = naexch + nbexch;

    if (nexch==0) {
      assert(make_pair(bra.alpha,bra.beta)==make_pair(ket.alpha,ket.beta));
      for (int r = 0; r < rnorb_; ++r) {
        out += (*Qaa)(bra.state, ket.state, r, r) * static_cast<double>(bra.alpha[r]);
        out += (*Qbb)(bra.state, ket.state, r, r) * static_cast<double>(bra.beta[r]);
      }
    }
    else if (nexch==2) {
      shared_ptr<const btas::Tensor3<double>> Q = (naexch==2 ? Qaa : Qbb);
      const int r = bit_to_numbers(naexch==2 ? (aexch & bra.alpha) : (bexch & bra.beta)).front();
      const int s = bit_to_numbers(naexch==2 ? (aexch & ket.alpha) : (bexch & ket.beta)).front();
      const double phase = static_cast<double>(bagel::sign((naexch==2 ? bra.alpha : bra.beta), r, s));
      out += phase * (*Q)(bra.state, ket.state, r, s);
    }
  }

  return out;
}

ProductCIHamTask::ProductCIHamTask(vector<PCI::Basis>* b, shared_ptr<const BlockOperators> blockops, shared_ptr<const DimerJop> jop, shared_ptr<const Matrix> mo1e, const size_t c1, double* d1, const size_t c2, double* d2) :
  ProductCITask<ProductCIHamTask>(b, c1, d1, c2, d2), blockops_(blockops), jop_(jop), mo1e_(mo1e), rnorb_(jop->monomer_jop<0>()->nocc()) {}

ProductCIHamiltonian::ProductCIHamiltonian(vector<PCI::Basis>& b, shared_ptr<const BlockOperators> blockops, shared_ptr<const DimerJop> jop) : Matrix(b.size(), b.size()), basis_(b), blockops_(blockops), jop_(jop) {
  const size_t size = basis_.size();

  shared_ptr<Matrix> mo1e = jop_->mo1e()->matrix();
  if (!jop_->hz()) {
    const int nocc = mo1e->mdim();
    for (int i = 0; i < nocc; ++i) {
      for (int j = 0; j < nocc; ++j) {
        for (int k = 0; k < nocc; ++k) {
          int a, b, c, d;
          tie(a,b) = minmax(i,k);
          tie(c,d) = minmax(k,j);
          mo1e->element(j,i) += 0.5*jop_->mo2e(a,b,c,d);
        }
      }
    }
  }

  const size_t ntasks = (size - 1)/2 + 1;
  TaskQueue<ProductCIHamTask> tasks(ntasks);

  size_t start = 0;
  size_t end = size - 1;

  for (size_t i = 0; i < ntasks; ++i, ++start, --end)
    tasks.emplace_back(&basis_, blockops_, jop_, mo1e, start, this->element_ptr(start, start), end, this->element_ptr(end, end));

  tasks.compute();

  this->fill_upper();
}
