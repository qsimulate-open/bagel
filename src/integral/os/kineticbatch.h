//
// BAGEL - Parallel electron correlation program.
// Filename: kineticbatch.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#ifndef __SRC_INTEGRAL_OS_KINETICBATCH_H
#define __SRC_INTEGRAL_OS_KINETICBATCH_H

#include <src/integral/os/osint.h>
#include <memory>

namespace bagel {

class KineticBatch : public OSInt {
  protected:
    void perform_VRR(double*) override;

  public:
    KineticBatch(const std::array<std::shared_ptr<const Shell>,2>&);
    ~KineticBatch();

    void compute() override;
};

}

#endif
