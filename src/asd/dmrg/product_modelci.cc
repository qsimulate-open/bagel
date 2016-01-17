//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: product_modelci.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <iostream>
#include <algorithm>

#include <src/asd/dmrg/product_modelci.h>
#include <src/asd/dmrg/product_citask.h>

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
  } else if (nexch == 2) {
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
  } else if (nexch == 4) {
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

double ProductCIHamTask::matrix_element_impl(PCI::Basis bra, PCI::Basis ket) {
  assert((bra.nelea+bra.neleb+bra.alpha.count()+bra.beta.count()) == (ket.nelea+ket.neleb+ket.alpha.count()+ket.beta.count()));

  const int ntransa = bra.alpha.count() - ket.alpha.count();
  const int ntransb = bra.beta.count() - ket.beta.count();

  const bitset<nbit__> aexch = bra.alpha ^ ket.alpha;
  const bitset<nbit__> bexch = bra.beta ^ ket.beta;

  const int naexch = aexch.count();
  const int nbexch = bexch.count();
  const int nexch = naexch + nbexch;

  double out = 0.0;

  // quick return for non-interacting determinants
  if (nexch > 4 || (std::abs(ntransa)+std::abs(ntransb)) > 2)
    return out;

  if (ntransa==0 && ntransb==0) {// "Diagonal contribution"
    if (bra.state==ket.state)
      out += compute_pure_ras(bra.alpha, bra.beta, ket.alpha, ket.beta);

    if (make_pair(bra.alpha,bra.beta)==make_pair(ket.alpha,ket.beta)) // |L>==|L'> --> <phi| H |phi'>
      out += blockops_->ham(ket.key(), bra.state, ket.state);

    if (nexch==0) {
      assert(make_pair(bra.alpha,bra.beta)==make_pair(ket.alpha,ket.beta));
      for (int r = 0; r < rnorb_; ++r) {
        out += blockops_->Q_aa(ket.key(), bra.state, ket.state, r, r) * static_cast<double>(bra.alpha[r]);
        out += blockops_->Q_bb(ket.key(), bra.state, ket.state, r, r) * static_cast<double>(bra.beta[r]);
      }
    } else if (nexch==2) {
      const int r = bit_to_numbers(naexch==2 ? (aexch & bra.alpha) : (bexch & bra.beta)).front();
      const int s = bit_to_numbers(naexch==2 ? (aexch & ket.alpha) : (bexch & ket.beta)).front();
      const double phase = sign((naexch==2 ? bra.alpha : bra.beta), r, s);
      if (naexch==2)
        out += phase * blockops_->Q_aa(ket.key(), bra.state, ket.state, r, s);
      else
        out += phase * blockops_->Q_bb(ket.key(), bra.state, ket.state, r, s);
    }
  } else if (ntransa*ntransb==-1) { // single spinflip
    if (naexch==1 && nbexch==1) {
      int r = bit_to_numbers(aexch).front();
      int s = bit_to_numbers(bexch).front();
      if (ntransa==-1)
        swap(bra,ket);
      const double phase = sign(ket.beta, s) * sign(ket.alpha, r) * static_cast<int>(1 - ((ket.alpha.count()%2) << 1));
      out += phase * blockops_->Q_ab(ket.key(), bra.state, ket.state, r, s);
    }
  } else if (ntransa*ntransb==1) { // alpha-beta transfer
    if (naexch==1 && nbexch==1) {
      int r = bit_to_numbers(aexch).front();
      int s = bit_to_numbers(bexch).front();
      if (ntransa==-1)
        swap(bra,ket);

      const double phase = sign(ket.beta, s) * sign(ket.alpha, r) * static_cast<int>(1 - ((ket.alpha.count()%2) << 1));
      out += phase * blockops_->P_ab(ket.key(), bra.state, ket.state, r, s);
    }
  } else if ((ntransa+ntransb)*(ntransa+ntransb)==1) { // ET
    if (nexch==1) {
      if ((ntransa+ntransb) > 0)
        swap(bra,ket);

      const int p = (abs(ntransa)==1 ? bit_to_numbers(aexch).front() : bit_to_numbers(bexch).front());
      const int etphase = 1 - (((bra.alpha.count()+bra.beta.count())%2) << 1);
      const double phase = etphase * (abs(ntransa)==1 ? sign(ket.alpha, p) : sign(ket.beta, p) * static_cast<int>(1 - ((ket.alpha.count()%2) << 1)));

      // 1e contribution
      out += phase * (abs(ntransa)==1 ? blockops_->S_a(ket.key(), bra.state, ket.state, p) : blockops_->S_b(ket.key(), bra.state, ket.state, p));

      bitset<nbit__> Jbit = abs(ntransa)==1 ? ket.beta : ket.alpha;
      bitset<nbit__> Kbit = abs(ntransa)==1 ? bra.alpha : bra.beta;
      if (abs(ntransa)==1)
        for (int r = 0; r < rnorb_; ++r)
          out += phase * (blockops_->D_a(ket.key(), bra.state, ket.state, r, r, p) * static_cast<double>(Jbit[r] + Kbit[r])
                           - blockops_->D_a(ket.key(), bra.state, ket.state, r, p, r) * static_cast<double>(Kbit[r]));
      else
        for (int r = 0; r < rnorb_; ++r)
          out += phase * (blockops_->D_b(ket.key(), bra.state, ket.state, r, r, p) * static_cast<double>(Jbit[r] + Kbit[r])
                           - blockops_->D_b(ket.key(), bra.state, ket.state, r, p, r) * static_cast<double>(Kbit[r]));
    } else if (nexch==3) {
      if ((ntransa+ntransb) > 0)
        swap(bra,ket);

      const int etphase = 1 - (((bra.alpha.count()+bra.beta.count())%2) << 1);
      if (naexch*nbexch==0) { // a^+a a or b^+b b
        assert(naexch==3 || nbexch==3);
        vector<int> annihilated = bit_to_numbers(naexch==3 ? aexch&ket.alpha : bexch&ket.beta);
        const int c = annihilated.back();
        const int b = annihilated.front();
        const int a = bit_to_numbers(naexch==3 ? aexch&bra.alpha : bexch&bra.beta).front();
        const double phase = etphase * (naexch==3 ? sign(bra.alpha, a, b) * sign(bra.alpha ^ aexch, c)
                                                  : sign(bra.beta, a, b) * sign(bra.beta ^ bexch, c) * static_cast<int>(1 - ((bra.alpha.count()%2) << 1)));
        if (naexch==3)
          out += phase * (blockops_->D_a(ket.key(), bra.state, ket.state, a, b, c) - blockops_->D_a(ket.key(), bra.state, ket.state, a, c, b));
        else
          out += phase * (blockops_->D_b(ket.key(), bra.state, ket.state, a, b, c) - blockops_->D_b(ket.key(), bra.state, ket.state, a, c, b));
      }
      else { // b^+b a or a^+a b
        if (make_pair(naexch, nbexch)==make_pair(2,1)) {
          const int c = bit_to_numbers(bexch).front();
          const int a = bit_to_numbers(aexch&ket.alpha).front();
          const int b = bit_to_numbers(aexch&bra.alpha).front();
          const double phase = etphase * sign(ket.alpha, a, b) * sign(ket.beta, c) * static_cast<int>(1 - ((bra.alpha.count()%2) << 1));
          out += phase * blockops_->D_b(ket.key(), bra.state, ket.state, a, b, c);
        }
        else if (make_pair(naexch, nbexch)==make_pair(1,2)) {
          const int c = bit_to_numbers(aexch).front();
          const int a = bit_to_numbers(bexch&ket.beta).front();
          const int b = bit_to_numbers(bexch&bra.beta).front();
          const double phase = etphase * sign(ket.beta, a, b) * sign(ket.alpha, c);
          out += phase * blockops_->D_a(ket.key(), bra.state, ket.state, a, b, c);
        }
      }
    }
  } else if ((ntransa+ntransb)*(ntransa+ntransb)==4) { // double ET
    if (nexch==2) {
      if ((ntransa+ntransb) < 0)
        swap(bra,ket);

      vector<int> exch_bits = bit_to_numbers(abs(ntransa)==2 ? aexch : bexch);
      const int p = exch_bits.front();
      const int q = exch_bits.back();
      const double phase = static_cast<int>(sign((abs(ntransa)==2 ? ket.alpha : ket.beta), q, p));
      if (abs(ntransa)==2)
        out += phase * (blockops_->P_aa(ket.key(), bra.state, ket.state, p, q) - blockops_->P_aa(ket.key(), bra.state, ket.state, q, p));
      else
        out += phase * (blockops_->P_bb(ket.key(), bra.state, ket.state, p, q) - blockops_->P_bb(ket.key(), bra.state, ket.state, q, p));
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
