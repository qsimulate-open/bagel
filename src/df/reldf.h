//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reldf.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
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


#ifndef __SRC_DF_RELDF_H
#define __SRC_DF_RELDF_H

#include <src/df/df.h>
#include <src/df/reldfbase.h>
#include <src/df/relcdmatrix.h>

namespace bagel {

class RelDFHalf;
class RelCDMatrix;

class RelDF : public RelDFBase, public std::enable_shared_from_this<RelDF> {
  protected:
    std::vector<int> alpha_;
    std::array<std::shared_ptr<const DFDist>,2> dfdata_;
    bool swap_;

    void set_basis() {
      auto ab0 = cartesian_.first  == Comp::L ? std::array<int,2>{{ Basis::LP, Basis::LM }}
                                              : std::array<int,2>{{ Basis::SP, Basis::SM }};
      auto ab1 = cartesian_.second == Comp::L ? std::array<int,2>{{ Basis::LP, Basis::LM }}
                                              : std::array<int,2>{{ Basis::SP, Basis::SM }};
      for (auto& i : ab0)
        for (auto& j : ab1)
          for (auto& a : alpha_) {
            auto tmp = std::make_shared<const SpinorInfo>(std::make_pair(i,j), a);
            if (tmp->nonzero(cartesian_)) basis_.push_back(tmp);
          }
    }

  public:
    RelDF(std::shared_ptr<const DFDist>, std::pair<int, int>, const std::vector<int>);
    RelDF(const RelDF&) = delete;
    RelDF(const RelDF&, bool);
    RelDF() = delete;

    std::array<std::shared_ptr<const DFDist>, 2> get_data() const { return dfdata_; }
    std::shared_ptr<const DFDist> get_real() const { return dfdata_[0]; }
    std::shared_ptr<const DFDist> get_imag() const { return dfdata_[1]; }

    bool not_diagonal() const { return cartesian_.first != cartesian_.second; }
    bool swapped() const { return swap_; }
    std::shared_ptr<const RelDF> swap() const;

    std::vector<std::shared_ptr<RelDFHalf>>
        compute_half_transform(std::array<std::shared_ptr<const Matrix>,4> r,
                               std::array<std::shared_ptr<const Matrix>,4> i) const;

    std::vector<std::shared_ptr<ZMatrix>> compute_Jop(std::list<std::shared_ptr<const RelCDMatrix>>& cd) const;
};

}

#endif
