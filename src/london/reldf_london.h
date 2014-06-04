//
// BAGEL - Parallel electron correlation program.
// Filename: reldf_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_LONDON_RELDF_LONDON_H
#define __SRC_LONDON_RELDF_LONDON_H

#include <src/rel/alpha.h>
#include <src/math/zmatrix.h>
#include <src/rel/reldfbase.h>
#include <src/london/cdmatrix_london.h>
#include <src/london/reldfhalf_london.h>

namespace bagel {

class RelDFHalf_London;
class CDMatrix_London;

class RelDF_London : public RelDFBase, public std::enable_shared_from_this<RelDF_London> {
  protected:
    std::vector<int> alpha_;
    std::shared_ptr<const ComplexDFDist> dfdata_;
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
    RelDF_London(std::shared_ptr<const ComplexDFDist>, std::pair<int, int>, const std::vector<int>);
    RelDF_London(const RelDF_London&) = delete;
    RelDF_London(const RelDF_London&, bool);
    RelDF_London() = delete;

    std::shared_ptr<const ComplexDFDist> df() const { return dfdata_; }
    bool not_diagonal() const { return cartesian_.first != cartesian_.second; }
    bool swapped() const { return swap_; }
    std::shared_ptr<const RelDF_London> swap() const;

    std::vector<std::shared_ptr<RelDFHalf_London>>
        compute_half_transform(std::array<std::shared_ptr<const Matrix>,4> r,
                               std::array<std::shared_ptr<const Matrix>,4> i) const;

    std::vector<std::shared_ptr<ZMatrix>> compute_Jop(std::list<std::shared_ptr<const CDMatrix_London>>& cd) const;
};

}

#endif
