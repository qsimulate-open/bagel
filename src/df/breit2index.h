//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: breit2index.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
// Maintainer: Shiozaki group
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


#ifndef __SRC_DF_BREIT2INDEX_H
#define __SRC_DF_BREIT2INDEX_H

#include <src/util/math/zmatrix.h>
#include <src/wfn/geometry.h>
#include <src/mat1e/rel/breitint.h>

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

