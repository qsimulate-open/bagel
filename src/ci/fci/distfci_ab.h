//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: distfci_ab.h
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

#ifndef __SRC_FCI_DISTFCI_AB_H
#define __SRC_FCI_DISTFCI_AB_H

#include <bitset>
#include <memory>
#include <src/util/f77.h>
#include <src/ci/fci/mofile.h>
#include <src/ci/fci/distcivec.h>

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

    std::vector<std::shared_ptr<RMATask<double>>> requests_;

  public:
    DistABTask(const std::bitset<nbit__> ast, std::shared_ptr<const Determinants> b, std::shared_ptr<const Determinants> i, std::shared_ptr<const MOFile> j,
               std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> s)
     : astring(ast), base_det(b), int_det(i), jop(j), sigma(s) {
      // first receive all the data (nele_a * lenb)
      const int norb_ = base_det->norb();
      const size_t lbs = base_det->lenb();

      buf = std::unique_ptr<double[]>(new double[lbs*(norb_-astring.count())]);

      for (int i = 0, k = 0; i != norb_; ++i) {
        if (!astring[i]) {
          std::bitset<nbit__> tmp = astring; tmp.set(i);
          auto j = cc->rma_rget(buf.get()+lbs*k++, base_det->lexical<0>(tmp));
          requests_.push_back(j);
        }
      }
    }

    void wait() {
      for (auto& i : requests_)
        i->wait();
    }

     std::list<std::shared_ptr<RMATask<double>>> compute() {
      const int norb_ = base_det->norb();
      const size_t lbs = base_det->lenb();
      const size_t lbt = int_det->lenb();

      const size_t ndima = norb_ - astring.count();

      btas::Tensor3<double> buf2(lbt, norb_, ndima);
      btas::Tensor3<double> buf3(lbt, norb_, ndima);
      btas::Tensor2<double> h(norb_*ndima, norb_*ndima);
      std::fill(buf2.begin(), buf2.end(), 0.0);

      // form Hamiltonian
      auto hiter = h.data();
      for (int i = 0; i != norb_; ++i) {
        if (astring[i]) continue;
        for (int j = 0; j != norb_; ++j) {
          for (int k = 0; k != norb_; ++k) {
            if (astring[k]) continue;
            for (int l = 0; l != norb_; ++l)
              *hiter++ = jop->mo2e(l+norb_*k,j+norb_*i);
          }
        }
      }

      for (int k = 0, i = 0; k != norb_; ++k) {
        if (!astring[k]) {
          for (int l = 0; l != norb_; ++l)
            for (auto& b : int_det->phiupb(l))
              buf2(b.source, l, i) += base_det->sign(astring, -1, k) * b.sign * buf[b.target+i*lbs];
          ++i;
        }
      }

      auto buf2v = btas::group(buf2, 1,3);
      auto buf3v = btas::group(buf3, 1,3);
      btas::contract(1.0, buf2v, {0,1}, h, {1,2}, 0.0, buf3v, {0,2});

      std::list<std::shared_ptr<RMATask<double>>> acctasks;
      for (int i = 0, k = 0; i < norb_; ++i) {
        if (astring[i]) continue;
        std::bitset<nbit__> atarget = astring; atarget.set(i);
        const double asign = base_det->sign(astring, -1, i);

        std::unique_ptr<double[]> bcolumn(new double[lbs]);
        std::fill_n(bcolumn.get(), lbs, 0.0);

        for (int j = 0; j < norb_; ++j) {
          for (auto& b : int_det->phiupb(j))
            bcolumn[b.target] += asign * b.sign * buf3(b.source, j, k);
        }
        std::shared_ptr<RMATask<double>> acctask = sigma->rma_radd(std::move(bcolumn), base_det->lexical<0>(atarget));
        acctasks.push_back(acctask);
        ++k;
      }
      return acctasks;
    }
};

}

#endif
