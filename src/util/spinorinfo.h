//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: spinorinfo.h
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

#ifndef __SRC_REL_SPINORINFO_H
#define __SRC_REL_SPINORINFO_H

#include <array>
#include <src/util/alpha.h>
#include <src/util/math/matop.h>

namespace bagel {

class SpinorInfo {
  protected:
    // X and Y (L+, L-, S+, S-)
    std::pair<int, int> basis_;
    // alpha component
    int alpha_comp_;
    int alpha_orig_;

  public:
    SpinorInfo(std::pair<int, int> moblock, const int alpha, const int alpha_orig = -1)
      : basis_(moblock), alpha_comp_(alpha), alpha_orig_(alpha_orig == -1 ? alpha : alpha_orig) {
    }

    std::complex<double> fac(const std::pair<int, int> cartesian) const {
      ZMatrix p1(4,1,true);
      ZMatrix p2(4,1,true);
      p1.element(basis_.first,  0) = 1.0;
      p2.element(basis_.second, 0) = 1.0;

      const Sigma s1(cartesian.first);
      const Sigma s2(cartesian.second);
      const Alpha a(alpha_orig_);

      ZMatrix z1(s1*p1);
      ZMatrix z2(s2*p2);
      return (z1 % a * z2).element(0,0);
    }
    bool nonzero(const std::pair<int, int> cartesian) const { return std::abs(fac(cartesian)) > 1.0e-20; }

    bool operator==(const SpinorInfo& o) const { return basis_.first == o.basis_.first && basis_.second == o.basis_.second && alpha_comp_ == o.alpha_comp_; }
    bool operator!=(const SpinorInfo& o) const { return !(*this == o); }

    std::shared_ptr<const SpinorInfo> swap() const {
      auto out = std::make_shared<SpinorInfo>(*this);
      std::swap(out->basis_.first, out->basis_.second);
      return out;
    }

    const std::pair<int, int>& basis() const { return basis_; }
    int basis(const int i) const {
      assert(i == 0 || i == 1);
      return i == 0 ? basis_.first : basis_.second;
    }

    int alpha_comp() const { return alpha_comp_; }
};

}

#endif
