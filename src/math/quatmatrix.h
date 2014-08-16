//
// BAGEL - Parallel electron correlation program.
// Filename: quatmatrix.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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
//
#ifndef __SRC_MATH_QUATMATRIX_H
#define __SRC_MATH_QUATMATRIX_H

#include <src/math/zmatrix.h>
#include <src/math/algo.h>

namespace bagel {

class QuatMatrix : public ZMatrix {
  public:
    QuatMatrix(const ZMatrix& o) : ZMatrix(o) { }
    QuatMatrix(ZMatrix&& o) : ZMatrix(std::move(o)) { }
    // TODO : implement symmetry checking function, and constructor that can take "00" and "01" matrices and build the rest of the matrix

    void diagonalize(VecView eig) override {
      assert(ndim() == mdim());
      // TODO parallelize
      zquatev_(ndim(), data(), eig.data()); 
      synchronize();
    }

    void synchronize() {
      mpi__->broadcast(data(), size(), 0);
    }

};

}

#endif
