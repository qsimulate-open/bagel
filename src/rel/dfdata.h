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

namespace bagel {

class DFData {
  protected:
    std::shared_ptr<const DFDist> dfdata_;
    std::pair<int, int> coord_;
    std::pair<int, int> basis_; 
    bool swap_;

    DFData(const DFData&, bool , bool);

  public:
    DFData(std::shared_ptr<const DFDist>, std::pair<int, int>);
    DFData(const DFData&) = delete;
    DFData() = delete;

    std::shared_ptr<const DFDist> df() const { return dfdata_; }
    std::pair<int, int> coord() const { return coord_; }
    std::pair<int, int> basis() const { return basis_; }
    bool cross() const { return coord_.first != coord_.second; }
    bool swapped() const { return swap_; }
    double cross_coeff() const;
    const int coeff_index() const;
    std::shared_ptr<const DFData> opp();
    std::shared_ptr<const DFData> swap();
    std::shared_ptr<const DFData> opp_and_swap();

    // to make readable
    enum Comp { X = 0, Y = 1, Z = 2, L = 3 }; 
    enum Basis { a = 0, b = 1 }; 

};

}

#endif
