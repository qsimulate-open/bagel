//
// BAGEL - Parallel electron correlation program.
// Filename: cispin.cc
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

using SD = pair<bitset<nbit__>, bitset<nbit__>>;

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

        const bitset<nbit__> aexch = abra ^ aket;
        const bitset<nbit__> bexch = bbra ^ bket;
        const int naexch = aexch.count();
        const int nbexch = bexch.count();
        const int nexch = naexch + nbexch;

        const int norb = norb_;

        if (nexch == 0) {
          const double sz = 0.5*static_cast<double>(bra.first.count() - bra.second.count());
          out += sz*sz + sz + static_cast<double>(bket.count());
        }
        else if (nexch == 4) {
          vector<int> bra_indices = to_vector(bbra);
          vector<int> ket_indices = to_vector(bket);
          for (auto& ibra : bra_indices) {
            for (auto& jket : ket_indices) {
              bitset<nbit__> abrap = abra; abrap.set(ibra);
              bitset<nbit__> bbrap = bbra; bbrap.reset(ibra);
              bitset<nbit__> aketp = aket; aketp.set(jket);
              bitset<nbit__> bketp = bket; bketp.reset(jket);
              if ( (abrap == aketp) && (bbrap == bketp) ) {
                const double phase = static_cast<double>(sign(bra.first, bra.second, ibra) * sign(ket.first, ket.second, jket));
                out += phase;
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

  this->fill_upper();
}
