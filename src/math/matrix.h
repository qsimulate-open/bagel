//
// BAGEL - Parallel electron correlation program.
// Filename: matrix.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_UTIL_MATRIX_H
#define __SRC_UTIL_MATRIX_H

#include <string>
#include <list>
#include <src/util/f77.h>
#include <src/math/matrix_base.h>
#include <src/math/distmatrix_base.h>

namespace bagel {

#ifdef HAVE_SCALAPACK
class DistMatrix;
#else
class Matrix;
using DistMatrix = Matrix;
#endif

class Matrix : public Matrix_base<double>, public std::enable_shared_from_this<Matrix> {
  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Matrix_base<double>>(*this);
    }

  public:
#ifdef HAVE_SCALAPACK
    Matrix(const int n, const int m, const bool localized = false);
#else
    Matrix(const int n, const int m, const bool localized = true);
#endif
    Matrix(const Matrix&);
    Matrix(const MatView&);
    Matrix(Matrix&&);
    Matrix() { }
    virtual ~Matrix() { }

    std::shared_ptr<Matrix> cut(const int nstart, const int nend) const { return get_submatrix(nstart, 0, nend-nstart, mdim()); }
    std::shared_ptr<Matrix> slice_copy(const int mstart, const int mend) const { return get_submatrix(0, mstart, ndim(), mend-mstart); }
    std::shared_ptr<Matrix> resize(const int n, const int m) const { return this->resize_impl<Matrix>(n, m); }
    std::shared_ptr<Matrix> merge(const std::shared_ptr<const Matrix> o) const { return this->merge_impl<Matrix>(o); }

    MatView slice(const int mstart, const int mend) {
      auto low = {0, mstart};
      auto up  = {ndim(), mend};
      assert(mstart >= 0 && mend <= mdim());
      return MatView(btas::make_rwview(this->range().slice(low, up), this->storage()), localized_);
    }

    const MatView slice(const int mstart, const int mend) const {
      auto low = {0, mstart};
      auto up  = {ndim(), mend};
      assert(mstart >= 0 && mend <= mdim());
      return MatView(btas::make_rwview(this->range().slice(low, up), this->storage()), localized_);
    }

    // antisymmetrize
    void antisymmetrize();

    // diagonalize this matrix (overwritten by a coefficient matrix)
    void diagonalize(VecView vec) override;
    std::shared_ptr<Matrix> diagonalize_blocks(VectorB& eig, std::vector<int> blocks) { return diagonalize_blocks_impl<Matrix>(eig, blocks); }
    std::tuple<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> svd(double* sing = nullptr);
    // compute S^-1. Assumes positive definite matrix
    void inverse();
    // compute S^-1 using symmetric form. Returns if some roots were discarded
    bool inverse_symmetric(const double thresh = 1.0e-8);
    // compute S^-1/2. If an eigenvalue of S is smaller than thresh, the root will be discarded. Returns if some roots were discarded
    bool inverse_half(const double thresh = 1.0e-8);
    std::shared_ptr<Matrix> tildex(const double thresh = 1.0e-8) const;
    // compute S^1/2. Same algorithm as above.
    void sqrt();

    using Matrix_base<double>::copy_block;
    using Matrix_base<double>::add_block;

    std::shared_ptr<Matrix> get_submatrix(const int nstart, const int mstart, const int ndim, const int mdim) const {
      return this->get_submatrix_impl<Matrix>(nstart, mstart, ndim, mdim);
    }

    Matrix& operator=(const Matrix& o) { Matrix_base<double>::operator=(o); return *this; }
    Matrix& operator=(Matrix&& o)      { Matrix_base<double>::operator=(o); return *this; }

    Matrix& operator/=(const Matrix&);
    Matrix operator/(const Matrix&) const;

    std::shared_ptr<Matrix> clone() const { return std::make_shared<Matrix>(ndim(), mdim(), localized_); }
    std::shared_ptr<Matrix> copy() const { return std::make_shared<Matrix>(*this); }

    // returns exp(*this)
    std::shared_ptr<Matrix> exp(const int deg = 6) const;
    // returns log(*this)
    std::shared_ptr<Matrix> log(const int deg = 6) const;
    // returns transpose(*this)
    std::shared_ptr<Matrix> transpose(const double a = 1.0) const;

    using Matrix_base<double>::ax_plus_y;
    using Matrix_base<double>::dot_product;
    void ax_plus_y(const double a, const Matrix& o) { this->ax_plus_y_impl(a, o); }
    double dot_product(const Matrix& o) const { return this->dot_product_impl(o); }

    double orthog(const std::list<std::shared_ptr<const Matrix>> o) { return this->orthog_impl(o); }
    void rotate(const int i, const int j, const double c, const double s) { drot_(ndim(), element_ptr(0,i), 1, element_ptr(0,j), 1, c, s); }
    void rotate(const int i, const int j, const double gamma) { rotate(i, j, cos(gamma), sin(gamma)); }
    void rotate(std::vector<std::tuple<int, int, double>>& rotations);

    // purify a (near unitary) matrix to be unitary
    void purify_unitary();
    void purify_idempotent(const Matrix& s);
    void purify_redrotation(const int nclosed, const int nact, const int nvirt);

    std::shared_ptr<Matrix> solve(std::shared_ptr<const Matrix> A, const int n) const;

#ifdef HAVE_SCALAPACK
    // return a shared pointer to this ifndef HAVE_SCALAPACK
    std::shared_ptr<DistMatrix> distmatrix() const;

    Matrix(const DistMatrix&);
#else
    std::shared_ptr<Matrix> distmatrix() { return shared_from_this(); }
    std::shared_ptr<const Matrix> distmatrix() const;
    std::shared_ptr<Matrix> matrix() { return shared_from_this(); }
    std::shared_ptr<const Matrix> matrix() const { return shared_from_this(); }
#endif
    std::shared_ptr<const Matrix> form_density_rhf(const int n, const int off = 0) const;
    std::shared_ptr<Matrix> swap_columns(const int i, const int iblock, const int j, const int jblock) const;
};


#ifdef HAVE_SCALAPACK
// Not to be confused with Matrix. DistMatrix is distributed and only supported when SCALAPACK is turned on. Limited functionality
class DistMatrix : public DistMatrix_base<double> {
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      std::shared_ptr<Matrix> mat = matrix();
      ar << mat;
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      std::shared_ptr<Matrix> mat;
      ar >> mat;
      DistMatrix tmp(*mat);
      ndim_ = tmp.ndim_;
      mdim_ = tmp.mdim_;
      desc_ = tmp.desc_;
      localsize_ = tmp.localsize_;
      local_ = std::unique_ptr<double[]>(new double[tmp.size()]);
      *this = tmp;
    }

  public:
    DistMatrix() { }
    DistMatrix(const int n, const int m);
    DistMatrix(const DistMatrix&);
    DistMatrix(DistMatrix&&);
    DistMatrix(const Matrix&);

    void diagonalize(VecView vec) override;

    DistMatrix operator*(const DistMatrix&) const;
    DistMatrix& operator*=(const DistMatrix&);
    DistMatrix operator%(const DistMatrix&) const; // caution
    DistMatrix operator^(const DistMatrix&) const; // caution
    DistMatrix operator+(const DistMatrix& o) const { DistMatrix out(*this); out += o; return out; }
    DistMatrix operator-(const DistMatrix& o) const { DistMatrix out(*this); out -= o; return out; }
    DistMatrix& operator+=(const DistMatrix& o) { ax_plus_y(1.0, o); return *this; }
    DistMatrix& operator-=(const DistMatrix& o) { ax_plus_y(-1.0, o); return *this; }
    DistMatrix& operator=(const DistMatrix& o);
    DistMatrix& operator=(DistMatrix&& o);

    std::shared_ptr<DistMatrix> clone() const { return std::make_shared<DistMatrix>(ndim_, mdim_); }

    using DistMatrix_base<double>::scale;
    using DistMatrix_base<double>::ax_plus_y;
    using DistMatrix_base<double>::dot_product;

    void ax_plus_y(const double a, const DistMatrix& o) { this->ax_plus_y_impl(a,o); }
    double dot_product(const DistMatrix& o) const { return this->dot_product_impl(o); }

    void rotate(std::vector<std::tuple<int, int, double>> rotations);
    void rotate(const int i, const int j, const double gamma);

    std::shared_ptr<Matrix> matrix() const;

    std::shared_ptr<const DistMatrix> form_density_rhf(const int n, const int off = 0) const;
};
#endif

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Matrix)
#ifdef HAVE_SCALAPACK
BOOST_CLASS_EXPORT_KEY(bagel::DistMatrix)
#endif

namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<Matrix, T>::value>::type> {
    typedef Matrix type;
  };
}

#endif
