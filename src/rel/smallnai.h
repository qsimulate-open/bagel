//
// BAGEL - Parallel electron correlation program.
// Filename: smallnai.h
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


#ifndef __SRC_REL_SMALLNAI_H
#define __SRC_REL_SMALLNAI_H

#include <src/util/zmatrix.h>
#include <src/scf/matrix1earray.h>

namespace bagel {

class SmallNAI : public Matrix1eArray<4> {
  protected:
    void init() override;

  public:
    SmallNAI(const std::shared_ptr<const Geometry>);

    void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int);

    void print(const std::string n = "") const override;

};

}

#endif
