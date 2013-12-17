//
// BAGEL - Parallel electron correlation program.
// Filename: modelci.cc
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

#include <iostream>
#include <algorithm>

#include <src/fci/modelci.h>
#include <src/fci/citask.h>

using namespace std;
using namespace bagel;

// ColumnTask computes two columns for the sake of load balancing
class CIHamTask : public CITask {
  protected:
    shared_ptr<const MOFile> jop_;
    shared_ptr<const CSymMatrix> mo1e_;

    double mo2e(int i, int j, int k, int l) {
      if (i>j) swap(i,j); if(k>l) swap(k,l);
      return jop_->mo2e(i,j,k,l);
    }

    double matrix_element(const SD& bra, const SD& ket) override {
      const bitset<nbit__> abra = bra.first;
      const bitset<nbit__> bbra = bra.second;
      const bitset<nbit__> aket = ket.first;
      const bitset<nbit__> bket = ket.second;

      const bitset<nbit__> aexch = abra ^ aket;
      const bitset<nbit__> bexch = bbra ^ bket;

      const int naexch = aexch.count();
      const int nbexch = bexch.count();
      const int nexch = naexch + nbexch;

      double out = 0.0;

      const int norb = norb_;

      if (nexch == 0) {
        // diagonal contribution
        vector<int> aoccs = to_vector(abra);
        vector<int> boccs = to_vector(bbra);

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
        vector<int> exch_ind = to_vector(naexch==2 ? aexch : bexch);
        const int i = exch_ind.front();
        const int j = exch_ind.back();
        const double signij = static_cast<double>( sign((naexch==2?abra:bbra), i, j) );

        // one-body
        out += signij * mo1e_->element(i,j);

        // two-body
        bitset<nbit__> acomm = abra & aket;
        bitset<nbit__> bcomm = bbra & bket;
        for (int k = 0; k < norb; ++k) {
          const double nj = static_cast<double>(acomm[k] + bcomm[k]);
          out += nj * signij * jop_->mo2e(i, j, k, k);
        }

        for (int& k : to_vector((naexch==2?acomm:bcomm)))
          out -= signij * mo2e(i, k, j, k);
      }
      else if (nexch == 4) {
        // double exchange (two-body only)
        if (naexch == nbexch) {
          vector<int> alphas = to_vector(aexch);
          vector<int> betas  = to_vector(bexch);
          const int i = alphas.front();
          const int j = alphas.back();
          const int k = betas.front();
          const int l = betas.back();

          const int signij = sign(abra, i, j);
          const int signkl = sign(bbra, k, l);
          out += signij * signkl * jop_->mo2e(i,j,k,l);
        }
        else {
          bitset<nbit__> exch = (naexch==4 ? aexch : bexch);
          bitset<nbit__> eket = (naexch==4 ? aket : bket);
          bitset<nbit__> ebra = (naexch==4 ? abra : bbra);
          vector<int> ann_list = to_vector(exch & eket);
          vector<int> cre_list = to_vector(exch & ebra);
          const int i = ann_list.front();
          const int j = ann_list.back();
          const int k = cre_list.front();
          const int l = cre_list.back();
          bitset<nbit__> tmp = eket;
          tmp.reset(i); tmp.reset(j);
          const double phase = static_cast<double>(sign(eket, i, j) * sign(tmp, k, l));
          out += phase * (mo2e(i,k,j,l) - mo2e(i,l,k,j));
        }
      }

      return out;
    }

  public:
    CIHamTask(vector<SD>* b, shared_ptr<const MOFile>& jop, shared_ptr<const CSymMatrix>& mo1e, const size_t c1, double* d1, const size_t c2, double* d2) :
      CITask(b, jop->nocc(), c1, d1, c2, d2), jop_(jop), mo1e_(mo1e) {}

};

CIHamiltonian::CIHamiltonian(vector<SD> b, shared_ptr<const MOFile> jop) : Matrix(b.size(), b.size()), basis_(b), jop_(jop) {
  const size_t size = basis_.size();

  shared_ptr<const CSymMatrix> mo1e;
  if (jop_->hz()) mo1e = jop_->mo1e();
  else {
    shared_ptr<CSymMatrix> tmp = jop_->mo1e()->copy();
    const int nocc = tmp->nocc();
    for (int i = 0; i < nocc; ++i) {
      for (int j = 0; j <= i; ++j) {
        for (int k = 0; k < nocc; ++k) {
          int a, b, c, d;
          tie(a,b) = minmax(i,k);
          tie(c,d) = minmax(k,j);
          tmp->element(j,i) += 0.5*jop_->mo2e(a,b,c,d);
        }
      }
    }
    mo1e = tmp;
  }

  const size_t ntasks = (size - 1)/2 + 1;
  TaskQueue<CIHamTask> tasks(ntasks);

  size_t start = 0;
  size_t end = size - 1;

  for (size_t i = 0; i < ntasks; ++i, ++start, --end)
    tasks.emplace_back(&basis_, jop_, mo1e, start, this->element_ptr(start, start), end, this->element_ptr(end, end));

  tasks.compute();

  this->fill_upper();
}

// ColumnTask computes two columns for the sake of load balancing
class CISpinTask : public CITask {
  protected:
    int sign(const bitset<nbit__>& bit1, const bitset<nbit__>& bit2, const int& i) const {
      const bitset<nbit__> mask((1ull << i) - 1);
      const int n = (mask & bit1).count() + (mask & bit2).count();
      return (1 - 2*(n%2));
    }

    double matrix_element(const SD& bra, const SD& ket) override {
      const bitset<nbit__> openbra = bra.first ^ bra.second;
      const bitset<nbit__> openket = ket.first ^ ket.second;
      const bitset<nbit__> closedbra = bra.first & bra.second;
      const bitset<nbit__> closedket = ket.first & ket.second;

      double out = 0.0;
      if (openbra == openket && closedbra == closedket) {
        const bitset<nbit__> abra = bra.first & openbra;
        const bitset<nbit__> bbra = bra.second & openbra;
        const bitset<nbit__> aket = ket.first & openket;
        const bitset<nbit__> bket = ket.second & openket;

        const int nexch = (abra ^ aket).count() + (bbra ^ bket).count();

        const int norb = norb_;

        if (nexch == 0) {
          const double sz = 0.5*static_cast<double>(bra.first.count() - bra.second.count());
          out += sz*sz + sz + static_cast<double>(bket.count());
        }
        else if (nexch == 4) {
          vector<int> bra_indices = to_vector(bbra);
          vector<int> ket_indices = to_vector(bket);
          for (auto& ibra : bra_indices) {
            bitset<nbit__> abrap = abra; abrap.set(ibra);
            bitset<nbit__> bbrap = bbra; bbrap.reset(ibra);
            for (auto& jket : ket_indices) {
              bitset<nbit__> aketp = aket; aketp.set(jket);
              bitset<nbit__> bketp = bket; bketp.reset(jket);
              if ( (abrap == aketp) && (bbrap == bketp) ) {
                out += static_cast<double>(sign(bra.first, bra.second, ibra) * sign(ket.first, ket.second, jket));
              }
            }
          }
        }
      }

      return out;
    }

  public:
    CISpinTask(vector<SD>* b, const int norb, const size_t c1, double* d1, const size_t c2, double* d2) :
      CITask(b, norb, c1, d1, c2, d2) {}

};

CISpin::CISpin(vector<SD>& b, const int norb) : Matrix(b.size(), b.size()), basis_(b) {
  const size_t size = basis_.size();

  const size_t ntasks = (size - 1)/2 + 1;
  TaskQueue<CISpinTask> tasks(ntasks);

  size_t start = 0;
  size_t end = size - 1;

  for (size_t i = 0; i < ntasks; ++i, ++start, --end)
    tasks.emplace_back(&basis_, norb, start, this->element_ptr(start, start), end, this->element_ptr(end, end));

  tasks.compute();

  fill_upper();
}
