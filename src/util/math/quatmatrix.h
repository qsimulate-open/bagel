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

#include <src/util/math/zmatrix.h>
#include <src/util/math/algo.h>

namespace bagel {

class QuatMatrix : public ZMatrix {
  protected:
    // Average out errors in time-reversal symmetry
    void t_symmetrize() {
      assert(mdim()%2 == 0 && ndim()%2 == 0);
      const size_t m = mdim()/2;
      const size_t n = ndim()/2;
      for (size_t i=0; i!=n; ++i) {
        for (size_t j=0; j!=m; ++j) {
          element(i, j) = 0.5*(element(i, j)+std::conj(element(n+i, m+j)));
          element(n+i, j) = 0.5*(element(n+i, j)-std::conj(element(i, m+j)));
          element(n+i, m+j) = std::conj(element(i, j));
          element(i, m+j) = -std::conj(element(n+i, j));
        }
      }
    }

    // Check that the matrix is symmetric under time-reversal
    bool is_t_symmetric(const double thresh = 1.0e-6) const {
      assert(mdim()%2 == 0 && ndim()%2 == 0);
      const int m = mdim()/2;
      const int n = ndim()/2;

      std::shared_ptr<ZMatrix> u = get_submatrix(0, 0, n, m);
      std::shared_ptr<ZMatrix> v = get_submatrix(n, 0, n, m);

      u->ax_plus_y(-1.0, *get_submatrix(n, m, n, m)->get_conjg());
      v->ax_plus_y( 1.0, *get_submatrix(0, m, n, m)->get_conjg());

      const double err = u->rms() + v->rms();
      return err < thresh;
    }

  public:
    QuatMatrix(const ZMatrix& o) : ZMatrix(o) { assert(is_t_symmetric()); }
    QuatMatrix(ZMatrix&& o) : ZMatrix(std::move(o)) { assert(is_t_symmetric()); }
    // TODO : implement constructor that can take "00" and "01" matrices and build the rest of the matrix

    void diagonalize(VecView eig) override {
      assert(ndim() == mdim());
      assert(eig.size() >= ndim());
      // assert that matrix is hermitian to ensure real eigenvalues
      assert(is_hermitian(1.0e-10));

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
