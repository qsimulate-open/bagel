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


#ifndef __SRC_REL_BREIT2INDEX_H
#define __SRC_REL_BREIT2INDEX_H

#include <src/math/zmatrix.h>
#include <src/wfn/geometry.h>
#include <src/rel/breitint.h>

namespace bagel {

/* class for J^{-1/2} B J^{-1/2} */

class Breit2Index {
  protected:
    std::pair<const int, const int> index_;
    std::shared_ptr<const Matrix> data_;

  public:
    Breit2Index(std::pair<const int, const int>, std::shared_ptr<const Matrix> breit, std::shared_ptr<const Matrix> data2);
    Breit2Index(std::pair<const int, const int> index, std::shared_ptr<const Matrix> data) : index_(index), data_(data) { }

    std::shared_ptr<const Matrix> data() const { return data_; }
    const std::pair<const int, const int>& index() const { return index_; }

    /// returning the same quantity with swapped indices
    std::shared_ptr<Breit2Index> cross() const;

    void print() const;

};

}

#endif

