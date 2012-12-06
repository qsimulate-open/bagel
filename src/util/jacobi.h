//
// BAGEL - Parallel electron correlation program.
// Filename: jacobi.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __BAGEL_UTIL_JACOBI_H
#define __BAGEL_UTIL_JACOBI_H

#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <list>
#include <src/scf/geometry.h>
#include <src/util/matrix.h>
#include <src/util/constants.h>

namespace bagel {

class Jacobi_base {
  protected:
    std::shared_ptr<Matrix> Q_; // Stores the eigenvectors

    const int nbasis_;

  public:
    Jacobi_base(std::shared_ptr<Matrix> Q) : Q_(Q), nbasis_(Q->ndim()) {};

    virtual void rotate(const int k, const int l) = 0;
    void sweep();
};

/************************************************************
* JacobiDiag provides the routines that would be            *
* used to diagonalize the matrix A_ and stores the eigen-   *
* vectors into Q_(if provided).                             *
************************************************************/

class JacobiDiag : public Jacobi_base {
  protected:
    std::shared_ptr<Matrix> A_; // The matrix to be diagonalized

  public:
    JacobiDiag(std::shared_ptr<Matrix> A, std::shared_ptr<Matrix> Q) : Jacobi_base(Q), A_(A) {};

    void rotate(const int k, const int l) override;
};

}

#endif
