//
// BAGEL - Parallel electron correlation program.
// Filename: dfbreit.h
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


#ifndef __SRC_REL_DFBREIT_H
#define __SRC_REL_DFBREIT_H

#include <memory>
#include <string>
#include <map>
#include <src/df/df.h>
#include <src/wfn/reference.h>
#include <src/rel/alpha.h>
#include <src/util/zmatrix.h>
#include <src/rel/reldfbase.h>
#include <src/rel/dfhalfcomplex.h>

namespace bagel {

class DFBreit {
  protected:
    std::vector<std::shared_ptr<const Alpha>> alpha1_;
    std::vector<std::shared_ptr<const Alpha>> alpha2_;
    std::shared_ptr<const DFDist> dfbreit_;
    bool swap_;

    DFBreit(const DFBreit&, bool);

  public:
    DFBreit(std::shared_ptr<const DFDist>, const int, const int);
    DFBreit(const DFBreit&) = delete;
    DFBreit() = delete;

    std::shared_ptr<const DFDist> df() const { return dfbreit_; }

};

}

#endif
