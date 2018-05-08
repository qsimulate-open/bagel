//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multisite.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifndef __BAGEL_MULTISITE_MULTISITE_H
#define __BAGEL_MULTISITE_MULTISITE_H

#include <memory>
#include <vector>
#include <src/util/input/input.h>
#include <src/wfn/reference.h>

namespace bagel {

// MultiSite prepares coefficient for ASD-DMRG solver
class MultiSite {
  protected:
    std::shared_ptr<const Reference> hf_ref_;
    std::shared_ptr<const Reference> mref_;

    // system info
    std::vector<int> region_sizes_;     // number of atoms on each site; atoms should be ordered by sites in molecular geometry

    void localize(std::shared_ptr<const PTree> ldata);

  public:
    // constructor
    MultiSite(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref);

    // return functions
    std::shared_ptr<const Reference> mref() const { return mref_; }
};

}

#endif
