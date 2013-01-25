//
// BAGEL - Parallel electron correlation program.
// Filename: distfci_bb.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_FCI_DISTFCI_BB_H
#define __SRC_FCI_DISTFCI_BB_H

#include <bitset>
#include <src/util/constants.h>
#include <src/fci/determinants.h>
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
        // TODO buffer should be nelec size (not norb size)
        const size_t npack = norb_*(norb_-1)/2;
        std::unique_ptr<double[]> ints(new double[npack*la]);
        std::unique_ptr<double[]> ints2(new double[npack*la]);
        std::fill_n(ints.get(), npack*la, 0.0);

        // first gather elements with correct sign
        for (int i = 0, ij = 0; i != norb_; ++i) {
          if (b[i]) { ij += i; continue; }
          for (int j = 0; j < i; ++j, ++ij) {
            if(b[j]) continue;
            std::bitset<nbit__> bt = b; bt.set(i); bt.set(j);
            const double ij_phase = base_det->sign(bt,i,j);
            daxpy_(la, ij_phase, source+la*base_det->lexical<1>(bt), 1, ints.get()+la*ij, 1);
          }
        }

        // call dgemm
        dgemm_("N", "N", la, npack, npack, 1.0, ints.get(), la, hamil, npack, 0.0, ints2.get(), la);

        for (int k = 0, kl = 0; k != norb_; ++k) {
          if (b[k]) { kl += k; continue; }
          for (int l = 0; l < k; ++l, ++kl) {
            if (b[l]) continue;
            std::bitset<nbit__> bt = b; bt.set(k); bt.set(l);
            const double kl_phase = base_det->sign(bt,k,l);
            const size_t ib = base_det->lexical<1>(bt);
            std::lock_guard<std::mutex> lock((*mutex)[ib]);
            daxpy_(la, kl_phase, ints2.get()+la*kl, 1, target+la*ib, 1);
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
          daxpy_(la, ij_phase, source+la*base_det->lexical<1>(bt), 1, ints.get()+la*i, 1);
        }
        dgemm_("N", "N", la, norb_, norb_, 1.0, ints.get(), la, hamil, norb_, 0.0, ints2.get(), la); 
        for (int i = 0; i != norb_; ++i) {
          if (b[i]) continue;
          std::bitset<nbit__> bt = b; bt.set(i);
          const double ij_phase = base_det->sign(bt, -1, i);
          const size_t ib = base_det->lexical<1>(bt);
          std::lock_guard<std::mutex> lock((*mutex)[ib]);
          daxpy_(la, ij_phase, ints2.get()+la*i, 1, target+la*ib, 1);
        }
      } else {
        assert(false);
      }
    }
};

}
#endif
