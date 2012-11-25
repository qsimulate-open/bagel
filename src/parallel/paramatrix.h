//
// BAGEL - Parallel electron correlation program.
// Filename: paramatrix.h
// Copyright (C) 2012 Toru Shiozaki
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

#ifndef __SRC_PARALLEL_PARAMATRIX_H
#define __SRC_PARALLEL_PARAMATRIX_H

// parallel matrix, a derived class of Matrix
#include <src/parallel/scalapack.h>
#include <src/parallel/process.h>
#include <src/util/matrix.h>

namespace bagel {

class ParaMatrix : public Matrix {
  protected:

  public:
    ParaMatrix(const int n, const int m) : Matrix(n,m) {}
    ParaMatrix(const ParaMatrix& o) : Matrix(o) { }

    void allreduce();
    void broadcast(const int root = 0);

#ifdef HAVE_SCALAPACK
    void diagonalize(double* eig) override;
#endif
};

}

#endif
