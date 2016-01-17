//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: prop1etask.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __BAGEL_FCI_PROP1ETASKS_H
#define __BAGEL_FCI_PROP1ETASKS_H

#include <src/util/f77.h>

namespace bagel {

class Prop1eTask {
  protected:
    std::shared_ptr<const Civec> cc_;
    const std::bitset<nbit__> targetstring_;
    double* const target_;
    const double* const integrals_;

  public:
    Prop1eTask(std::shared_ptr<const Civec> cc, const std::bitset<nbit__>& targetstring, double* const target, const double* const integrals) :
      cc_(cc), targetstring_(targetstring), target_(target), integrals_(integrals) {}

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

          const double hc = integrals_[j + i*norb] * static_cast<double>(det->sign(sourcestring, i, j));
          blas::ax_plus_y_n(hc, cc_->element_ptr(0, det->lexical<0>(sourcestring)), lb, target_);
        }
      }
    }

};

}

#endif
