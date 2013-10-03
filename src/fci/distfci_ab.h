//
// BAGEL - Parallel electron correlation program.
// Filename: distfci_ab.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_FCI_DISTFCI_AB_H
#define __SRC_FCI_DISTFCI_AB_H

#include <bitset>
#include <memory>
#include <src/util/f77.h>
#include <src/fci/determinants.h>
#include <src/fci/mofile.h>
#include <src/fci/civec.h>

namespace bagel {

class DistABTask {
  protected:
    std::bitset<nbit__> astring;
    std::shared_ptr<const Determinants> base_det;
    std::shared_ptr<const Determinants> int_det;
    std::shared_ptr<const MOFile> jop;
    std::shared_ptr<const DistCivec> cc;
    std::shared_ptr<DistCivec> sigma;

    std::unique_ptr<double[]> buf;

    std::vector<int> requests_;

  public:
    DistABTask(const std::bitset<nbit__> ast, std::shared_ptr<const Determinants> b, std::shared_ptr<const Determinants> i, std::shared_ptr<const MOFile> j,
               std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> s)
     : astring(ast), base_det(b), int_det(i), jop(j), sigma(s) {
      // first receive all the data (nele_a * lenb)
      const int norb_ = base_det->norb();
      const size_t lbs = base_det->lenb();
      const size_t lbt = int_det->lenb();

      buf = std::unique_ptr<double[]>(new double[lbs*(norb_-astring.count())]);

      for (int i = 0, k = 0; i != norb_; ++i) {
        if (!astring[i]) {
          std::bitset<nbit__> tmp = astring; tmp.set(i);
          const int j = cc->get_bstring_buf(buf.get()+lbs*k++, base_det->lexical<0>(tmp));
          if (j >= 0) requests_.push_back(j);
        }
      }

    }

    bool test() {
      bool out = true;
      for (auto i = requests_.begin(); i != requests_.end(); )
        if (mpi__->test(*i)) {
          i = requests_.erase(i);
        } else {
          ++i;
          out = false;
        }
      return out;
    }

    void compute() {
      const int norb_ = base_det->norb();
      const size_t lbs = base_det->lenb();
      const size_t lbt = int_det->lenb();

      const size_t ndima = norb_ - astring.count();

      std::unique_ptr<double[]> buf2(new double[lbt*ndima*norb_]);
      std::unique_ptr<double[]> buf3(new double[lbt*ndima*norb_]);
      std::unique_ptr<double[]> h(new double[ndima*ndima*norb_*norb_]);
      std::fill_n(buf2.get(), lbt*ndima*norb_, 0.0);

      // form Hamiltonian
      for (int i = 0, ijkl = 0; i != norb_; ++i) {
        if (astring[i]) continue;
        for (int j = 0; j != norb_; ++j) {
          for (int k = 0; k != norb_; ++k) {
            if (astring[k]) continue;
            for (int l = 0; l != norb_; ++l)
              h[ijkl++] = jop->mo2e(l+norb_*k,j+norb_*i);
          }
        }
      }

      for (int k = 0, i = 0; k != norb_; ++k) {
        if (!astring[k]) {
          for (int l = 0; l != norb_; ++l)
            for (auto& b : int_det->phiupb(l))
              buf2[b.source+lbt*(l+norb_*i)] += base_det->sign(astring, -1, k) * b.sign * buf[b.target+i*lbs];
          ++i;
        }
      }

      const int ij = ndima*norb_;
      dgemm_("n", "n", lbt, ij, ij, 1.0, buf2, lbt, h, ij, 0.0, buf3, lbt);

      for (int i = 0, k = 0; i < norb_; ++i) {
        if (astring[i]) continue;
        std::bitset<nbit__> atarget = astring; atarget.set(i);
        const double asign = base_det->sign(astring, -1, i);

        std::unique_ptr<double[]> bcolumn(new double[lbs]);
        std::fill_n(bcolumn.get(), lbs, 0.0);

        for (int j = 0; j < norb_; ++j) {
          for (auto& b : int_det->phiupb(j))
            bcolumn[b.target] += asign * b.sign * buf3[b.source+lbt*(j+norb_*k)];
        }
        sigma->accumulate_bstring_buf(bcolumn, base_det->lexical<0>(atarget));
        ++k;
      }
      sigma->flush();
    }
};

}

#endif
