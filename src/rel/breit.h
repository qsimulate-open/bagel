//
// BAGEL - Parallel electron correlation program.
// Filename: breit.h
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


#ifndef __SRC_REL_BREIT_H
#define __SRC_REL_BREIT_H

#include <memory>
#include <array>
#include <src/util/zmatrix.h>
#include <src/util/matrix.h>
#include <src/wfn/geometry.h>
#include <src/rel/nmatrix1e.h>

namespace bagel {

class Breit : public NMatrix1e {
  protected:
    void init() override;

  public:
    Breit(const std::shared_ptr<const Geometry>);

    ~Breit() {};
  
    void computebatch(const std::array<std::shared_ptr<const Shell>,4>&, const int, const int);

    void print() const;

    // 6 blocks for breit
    constexpr static int nblocks() { return 6; }
};

}

#endif

