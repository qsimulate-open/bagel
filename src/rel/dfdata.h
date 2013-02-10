//
// BAGEL - Parallel electron correlation program.
// Filename: dfdata.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#ifndef __SRC_REL_DFDATA_H
#define __SRC_REL_DFDATA_H

#include <memory>
#include <string>
#include <map>
#include <src/df/df.h>
#include <src/wfn/reference.h>
#include <src/rel/alpha.h>
#include <src/util/zmatrix.h>
#include <src/rel/reldfbase.h>

namespace bagel {

class DFData : public RelDFBase {
  protected:
    std::shared_ptr<const DFDist> dfdata_;
    bool swap_;

    void set_basis() override {
      std::array<int, 2> ab = {{Basis::a, Basis::b}};
      for (auto& i : ab) {
        for (auto& j : ab) {
          std::shared_ptr<ABcases> tmp(new ABcases(std::make_pair(i,j), coord_, sigma1_, sigma2_, alpha_));
          if (tmp->nonzero()) basis_.push_back(tmp);
        }
      }
      assert(basis_.size() == 2);
    }

    DFData(const DFData&, bool);

  public:
    DFData(std::shared_ptr<const DFDist>, std::pair<int, int>, const int);
    DFData(const DFData&) = delete;
    DFData() = delete;

    std::shared_ptr<const DFDist> df() const { return dfdata_; }
    bool cross() const { return coord_.first != coord_.second; }
    bool swapped() const { return swap_; }
    std::shared_ptr<const DFData> swap() const;

//  const std::tuple<int, int, int, int> compute_index_Jop() const;

};

}

#endif
