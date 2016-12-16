//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: distfci_bb.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_FCI_DISTFCI_BB_H
#define __SRC_FCI_DISTFCI_BB_H

#include <bitset>
#include <src/util/constants.h>
#include <src/ci/fci/determinants.h>
#include <mutex>

namespace bagel {

class DistBBTask {
  protected:
    const size_t la;
    const double* source;
    double* target;
    const double* hamil;
    std::shared_ptr<const Determinants> base_det;
    const std::bitset<nbit__> b;
    std::vector<std::mutex>* mutex;

  public:
    DistBBTask(const size_t& l, const double* s, double* t, const double* h, std::shared_ptr<const Determinants> det, const std::bitset<nbit__>& bs, std::vector<std::mutex>* m)
      : la(l), source(s), target(t), hamil(h), base_det(det), b(bs), mutex(m) {}

    void compute() {
      const size_t norb_ = base_det->norb();
      if (b.count()+2 == base_det->neleb()) {
        const size_t nunoc = norb_ - b.count();
        const size_t npack = nunoc*(nunoc-1)/2;
        std::unique_ptr<double[]> ints(new double[npack*la]);
        std::unique_ptr<double[]> ints2(new double[npack*la]);
        std::fill_n(ints.get(), npack*la, 0.0);

        // form Hamiltonian
        std::unique_ptr<double[]> h(new double[npack*npack]);
        for (int i = 0, ijkl = 0, ij1 = 0; i != norb_; ++i) {
          if (b[i]) { ij1 += i; continue; }
          for (int j = 0; j < i; ++j, ++ij1) {
            if(b[j]) continue;
            for (int k = 0, kl1 = ij1*norb_*(norb_-1)/2; k != norb_; ++k) {
              if (b[k]) { kl1 += k; continue; }
              for (int l = 0; l < k; ++l, ++kl1)
                if (!b[l]) h[ijkl++] = hamil[kl1];
            }
          }
        }

        // first gather elements with correct sign
        for (int i = 0, ij = 0; i != norb_; ++i) {
          if (b[i]) continue;
          for (int j = 0; j < i; ++j) {
            if(b[j]) continue;
            std::bitset<nbit__> bt = b; bt.set(i); bt.set(j);
            const double ij_phase = base_det->sign(bt,i,j);
            blas::ax_plus_y_n(ij_phase, source+la*base_det->lexical<1>(bt), la, ints.get()+la*ij++);
          }
        }

        // call dgemm
        dgemm_("N", "N", la, npack, npack, 1.0, ints.get(), la, h.get(), npack, 0.0, ints2.get(), la);

        for (int k = 0, kl = 0; k != norb_; ++k) {
          if (b[k]) continue;
          for (int l = 0; l < k; ++l) {
            if (b[l]) continue;
            std::bitset<nbit__> bt = b; bt.set(k); bt.set(l);
            const double kl_phase = base_det->sign(bt,k,l);
            const size_t ib = base_det->lexical<1>(bt);
            std::lock_guard<std::mutex> lock((*mutex)[ib]);
            blas::ax_plus_y_n(kl_phase, ints2.get()+la*kl++, la, target+la*ib);
          }
        }
      } else if (b.count()+1 == base_det->neleb()) {
        std::unique_ptr<double[]> ints(new double[la*norb_]);
        std::unique_ptr<double[]> ints2(new double[la*norb_]);
        std::fill_n(ints.get(), la*norb_, 0.0);
        for (int i = 0; i != norb_; ++i) {
          if (b[i]) continue;
          std::bitset<nbit__> bt = b; bt.set(i);
          const double ij_phase = base_det->sign(bt, -1, i);
          blas::ax_plus_y_n(ij_phase, source+la*base_det->lexical<1>(bt), la, ints.get()+la*i);
        }
        dgemm_("N", "N", la, norb_, norb_, 1.0, ints.get(), la, hamil, norb_, 0.0, ints2.get(), la);
        for (int i = 0; i != norb_; ++i) {
          if (b[i]) continue;
          std::bitset<nbit__> bt = b; bt.set(i);
          const double ij_phase = base_det->sign(bt, -1, i);
          const size_t ib = base_det->lexical<1>(bt);
          std::lock_guard<std::mutex> lock((*mutex)[ib]);
          blas::ax_plus_y_n(ij_phase, ints2.get()+la*i, la, target+la*ib);
        }
      } else {
        assert(false);
      }
    }
};

}
#endif
