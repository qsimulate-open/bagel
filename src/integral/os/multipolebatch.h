//
// BAGEL - Parallel electron correlation program.
// Filename: multipolebatch.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_INTEGRAL_OS_MULTIPOLEBATCH_H
#define __SRC_INTEGRAL_OS_MULTIPOLEBATCH_H

#include <src/integral/os/multipolebatch_base.h>

namespace bagel {

class MultipoleBatch : public MultipoleBatch_base {
  protected:

  public:
    MultipoleBatch(const std::array<std::shared_ptr<const Shell>,2>& shells, const std::array<double, 3> centre,
                   const int lmax = ANG_HRR_END, std::shared_ptr<StackMem> = nullptr);
    ~MultipoleBatch() { }

    void compute() override;
};

}

#endif
