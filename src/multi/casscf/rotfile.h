//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rotfile.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __BAGEL_CASSCF_ROTFILE_H
#define __BAGEL_CASSCF_ROTFILE_H

#include <src/util/math/matrix.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

template<typename DataType>
class RotationMatrix {
  public:
    using data_type = DataType;
    using MatType  = typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type;
    using ViewType = typename std::conditional<std::is_same<DataType,double>::value, MatView, ZMatView>::type;

  protected:
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    const int size_;
    std::unique_ptr<DataType[]> data_;

  public:
    RotationMatrix(const int iclos, const int iact, const int ivirt);
    RotationMatrix(const RotationMatrix& o);
    RotationMatrix(std::shared_ptr<const RotationMatrix> o);
    RotationMatrix(std::shared_ptr<const Matrix_base<DataType>> o, const int iclos, const int iact, const int ivirt);

    std::shared_ptr<RotationMatrix<DataType>> clone() const;
    std::shared_ptr<RotationMatrix<DataType>> copy() const;

    // overloaded operators
    RotationMatrix<DataType> operator+(const RotationMatrix<DataType>& o) const;
    RotationMatrix<DataType> operator-(const RotationMatrix<DataType>& o) const;
    RotationMatrix<DataType>& operator+=(const RotationMatrix<DataType>& o);
    RotationMatrix<DataType>& operator-=(const RotationMatrix<DataType>& o);
    RotationMatrix<DataType>& operator*=(const DataType a);
    RotationMatrix<DataType>& operator/=(const RotationMatrix<DataType>& o);
    RotationMatrix<DataType>& operator*=(const RotationMatrix<DataType>& o);
    RotationMatrix<DataType> operator/(const RotationMatrix<DataType>& o) const;
    RotationMatrix<DataType> operator*(const RotationMatrix<DataType>& o) const;
    RotationMatrix<DataType>& operator=(const RotationMatrix<DataType>& o);

    // size of the file
    int size() const { return size_; }
    // zero out
    void zero() { fill(DataType(0.0)); }
    void fill(const DataType& a) { std::fill_n(data(), size_, a); }
    // returns dot product
    DataType dot_product(const RotationMatrix<DataType>& o) const { return blas::dot_product(data(), size_, o.data()); }
    DataType dot_product(const std::shared_ptr<const RotationMatrix<DataType>> o) const { return dot_product(*o); }
    // scale function
    void scale(const DataType& a) { std::for_each(data(), data()+size_, [&a](DataType& p) { p *= a; }); }
    // returns norm of the vector
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    // returns a root mean square
    double rms() const { return norm() / std::sqrt(static_cast<double>(size())); }
    // daxpy added to self
    void ax_plus_y(const DataType& a, const RotationMatrix& o) { blas::ax_plus_y_n(a, o.data(), size_, data()); }
    void ax_plus_y(const DataType& a, const std::shared_ptr<const RotationMatrix> o) { ax_plus_y(a, *o); }

    // orthogonalize to the liset of RotationMatrix's
    double orthog(std::list<std::shared_ptr<const RotationMatrix<DataType>>> c);
    double normalize();
    void synchronize();

    // return data_
    DataType* data() { return data_.get(); }
    const DataType* data() const { return data_.get(); }
    DataType& data(const size_t i) { return data_[i]; }
    const DataType& data(const size_t i) const { return data_[i]; }
    // return data_
    DataType* begin() { return data(); }
    // return data_
    DataType* end() { return data()+size_; }

    // closed-active block. closed runs first
    DataType* ptr_ca() { return data(); }
    DataType& ele_ca(const int ic, const int ia) { return data_[ic + ia*nclosed_]; }
    // active-virtual block. virtual runs first
    DataType* ptr_va() { return data() + nclosed_*nact_; }
    DataType& ele_va(const int iv, const int ia) { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    // closed-virtual block. virtual runs first
    DataType* ptr_vc() { return data() + (nclosed_+nvirt_)*nact_; }
    DataType& ele_vc(const int iv, const int ic) { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    // const references and pointers
    const DataType& ele_ca(const int ic, const int ia) const { return data_[ic + ia*nclosed_]; }
    const DataType& ele_va(const int iv, const int ia) const { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    const DataType& ele_vc(const int iv, const int ic) const { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    const DataType* ptr_ca() const { return data(); }
    const DataType* ptr_va() const { return data() + nclosed_*nact_; }
    const DataType* ptr_vc() const { return data() + (nclosed_+nvirt_)*nact_; }

    void ax_plus_y_ca(const DataType a, const ViewType mat);
    void ax_plus_y_va(const DataType a, const ViewType mat);
    void ax_plus_y_vc(const DataType a, const ViewType mat);

    std::shared_ptr<MatType> ca_mat() const;
    std::shared_ptr<MatType> va_mat() const;
    std::shared_ptr<MatType> vc_mat() const;

    // unpack to Matrix
    std::shared_ptr<MatType> unpack(const DataType a = 0.0) const;
    std::shared_ptr<MatType> unpack_sym(const DataType a = 0.0) const;

    // print matrix
    void print(const std::string in = "") const;
};

using RotFile = RotationMatrix<double>;
using ZRotFile = RotationMatrix<std::complex<double>>;

extern template class RotationMatrix<double>;
extern template class RotationMatrix<std::complex<double>>;

}

#endif
