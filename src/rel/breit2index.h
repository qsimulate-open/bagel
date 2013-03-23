//
// BAGEL - Parallel electron correlation program.
// Filename: breit2index.h
// Copyright (C) 2013 Matthew Kelley
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
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


#ifndef __SRC_REL_BREIT2INDEX_H
#define __SRC_REL_BREIT2INDEX_H

#include <memory>
#include <array>
#include <src/util/zmatrix.h>
#include <src/util/matrix.h>
#include <src/wfn/geometry.h>
#include <src/rel/breit.h>

namespace bagel {

class Breit2Index {
  protected:
    std::pair<const int, const int> index_;
    std::shared_ptr<const Matrix> k_term_;
    Breit2Index(std::pair<const int, const int>, std::shared_ptr<const Matrix> k);

  public:
    Breit2Index(std::pair<const int, const int>, std::shared_ptr<const Matrix> breit, std::shared_ptr<const Matrix> data2);

    std::shared_ptr<const Matrix> k_term() const { return k_term_; }
    const std::pair<const int, const int>& index() const { return index_; }

    std::shared_ptr<Breit2Index> cross() const;

    void print() const;

};

}

#endif

