//
// BAGEL - Parallel electron correlation program.
// Filename: HZtasks.h
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

#ifndef __BAGEL_FCI_HZTASKS_H
#define __BAGEL_FCI_HZTASKS_H

#include <src/util/f77.h>

namespace bagel {

template<typename DataType>
class HZTaskAA {
  protected:
    const std::shared_ptr<const Civector<DataType>> cc_;
    const std::bitset<nbit__> targetstring_;
    DataType* const target_;
    const DataType* const h1_;
    const DataType* const h2_;

  public:
    HZTaskAA(std::shared_ptr<const Civector<DataType>> cc, const std::bitset<nbit__> targetstring, DataType* const target,
             const DataType* const h1, const DataType* const h2) :
      cc_(cc), targetstring_(targetstring), target_(target), h1_(h1), h2_(h2) {}

    void compute() {
      std::shared_ptr<const Determinants> det = cc_->det();

      const int lb = cc_->lenb();
      const int norb = det->norb();

      // One-electron part
      for (int i = 0; i < norb; ++i) {
        if (!targetstring_[i]) continue;
        std::bitset<nbit__> ibs = targetstring_; ibs.reset(i);

        for (int j = 0; j < norb; ++j) {
          if (ibs[j]) continue;
          std::bitset <nbit__> sourcestring = ibs; sourcestring.set(j);

          const DataType hc = h1_[i+ j*norb] * static_cast<double>(det->sign(sourcestring, i, j));
          const DataType* const source = cc_->element_ptr(0, det->lexical<0>(sourcestring));
          blas::ax_plus_y_n(hc, source, lb, target_);
        }
      }

      // Two-electron part
      for (int i = 0; i != norb; ++i) {
        if (!targetstring_[i]) continue;
        for (int j = 0; j < i; ++j) {
          if(!targetstring_[j]) continue;
          const int ij_phase = det->sign(targetstring_,i,j);
          std::bitset<nbit__> string_ij = targetstring_;
          string_ij.reset(i); string_ij.reset(j);
          for (int l = 0; l != norb; ++l) {
            if (string_ij[l]) continue;
            for (int k = 0; k < l; ++k) {
              if (string_ij[k]) continue;
              const int kl_phase = det->sign(string_ij,l,k);
              const double phase = -static_cast<double>(ij_phase*kl_phase);
              std::bitset<nbit__> string_ijkl = string_ij;
              string_ijkl.set(k); string_ijkl.set(l);
              const DataType temp = phase * h2_[i+norb*(j+norb*(k+norb*l))]; // TODO fix this random access (but stick to this order for complex cases)
              const DataType* const source = cc_->element_ptr(0, det->lexical<0>(string_ijkl));
              std::transform(source, source+lb, target_, target_, [&temp](DataType p, DataType q){ return temp*p+q; });
            }
          }
        }
      }
    }

};


template<typename DataType>
class HZTaskAB1 {
  protected:
    std::shared_ptr<const Determinants> det_;
    const int lbs_;
    const DataType* const source_base_;
    DataType* const target_base_;
    const int k_;
    const int l_;


  public:
    HZTaskAB1(std::shared_ptr<const Determinants>& det, const int& lbs, const DataType* const source_base, DataType* const target_base,
      const int& k, const int& l) :
      det_(det), lbs_(lbs), source_base_(source_base), target_base_(target_base), k_(k), l_(l) {}

    void compute() {
      const int lbt = det_->lenb();

      for (auto& aiter : det_->phiupa(k_)) {
        DataType* target = target_base_ + aiter.source*lbt;
        const DataType* source = source_base_ + aiter.target*lbs_;
        for (auto& biter : det_->phiupb(l_)) {
          const double sign = aiter.sign * biter.sign;
          target[biter.source] += sign * source[biter.target];
        }
      }
    }

};

}

#endif
