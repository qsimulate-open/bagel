//
// BAGEL - Parallel electron correlation program.
// Filename: breitterm.h
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


#ifndef __SRC_REL_BREITTERM_H
#define __SRC_REL_BREITTERM_H

#include <memory>
#include <array>
#include <src/util/zmatrix.h>
#include <src/util/matrix.h>
#include <src/wfn/geometry.h>
#include <src/rel/dfdata.h>
#include <src/rel/breit.h>

namespace bagel {

class BreitTerm {
  protected:
    std::shared_ptr<const Breit> breit_;
    std::vector<std::shared_ptr<ZMatrix>> bt_;
    std::vector<std::pair<const int, const int>> index_;
    std::array<std::list<std::shared_ptr<ZMatrix>>, 6> data_;

  public:
    BreitTerm(std::shared_ptr<const Breit>, std::list<std::shared_ptr<DFData>>, std::list<std::shared_ptr<const ZMatrix>>, std::vector<int>);

    ~BreitTerm() {};

    std::list<std::shared_ptr<ZMatrix>> data(const int i) { return data_[i]; }
    std::array<std::list<std::shared_ptr<ZMatrix>>, 6> data() { return data_; }

};

}

#endif

