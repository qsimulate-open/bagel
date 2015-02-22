//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbital/smallrotfile.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
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

#ifndef __BAGEL_ASD_ACTIVEROTFILE_H
#define __BAGEL_ASD_ACTIVEROTFILE_H

#include <src/util/math/matrix.h>

namespace bagel {

template<typename DataType>
class ASD_ActiveRotationMatrix {
  public:
    using data_type = DataType;

  protected:
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    const int nactA_;
    const int nactB_;
    const int size_;
    std::unique_ptr<DataType[]> data_;

  public:
    ASD_ActiveRotationMatrix(const int iclos, const int iact, const int ivirt, const int iactA, const int iactB)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), nactA_(iactA), nactB_(iactB), size_(iactA*iactB), data_(new DataType[size_]) {
      zero();
    }
    ASD_ActiveRotationMatrix(const ASD_ActiveRotationMatrix& o) : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), nactA_(o.nactA_), nactB_(o.nactB_),  size_(o.size_), data_(new DataType[o.size_]) {
      *this = o;
    }
    ASD_ActiveRotationMatrix(std::shared_ptr<const ASD_ActiveRotationMatrix> o)
      : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), nactA_(o->nactA_), nactB_(o->nactB_), size_(o->size_), data_(new DataType[o->size_]) {
      *this = *o;
    }
    //Matrix to ASDRotationMatrix conversion
    ASD_ActiveRotationMatrix(std::shared_ptr<const Matrix_base<DataType>> o, const int iclos, const int iact, const int ivirt, const int iactA, const int iactB)
      : nclosed_(iclos), nact_(iact), nvirt_(ivirt), nactA_(iactA), nactB_(iactB), size_(iactA*iactB), data_(new DataType[size_]) {
    //const int nocc = nclosed_ + nact_;
      assert(iact == iactA+iactB);
      for (int i = 0; i != nactA_; ++i) {
        for (int j = 0; j != nactB_; ++j) {
          ele_aa(j, i) = o->element(j+nclosed_+nactA_,i+nclosed_);
        }
      }
    }

    std::shared_ptr<ASD_ActiveRotationMatrix<DataType>> clone() const {
      return std::make_shared<ASD_ActiveRotationMatrix<DataType>>(nclosed_, nact_, nvirt_, nactA_, nactB_);
    }
    std::shared_ptr<ASD_ActiveRotationMatrix<DataType>> copy() const {
      return std::make_shared<ASD_ActiveRotationMatrix<DataType>>(*this);
    }

    // overloaded operators
    ASD_ActiveRotationMatrix<DataType> operator+(const ASD_ActiveRotationMatrix<DataType>& o) const { ASD_ActiveRotationMatrix<DataType> out(*this); out.ax_plus_y(1.0, o); return out; }
    ASD_ActiveRotationMatrix<DataType> operator-(const ASD_ActiveRotationMatrix<DataType>& o) const { ASD_ActiveRotationMatrix<DataType> out(*this); out.ax_plus_y(-1.0, o); return out; }
    ASD_ActiveRotationMatrix<DataType>& operator+=(const ASD_ActiveRotationMatrix<DataType>& o) { ax_plus_y(1.0, o); return *this; }
    ASD_ActiveRotationMatrix<DataType>& operator-=(const ASD_ActiveRotationMatrix<DataType>& o) { ax_plus_y(-1.0, o); return *this; }
    ASD_ActiveRotationMatrix<DataType>& operator*=(const DataType a) { scale(a); return *this; }
    ASD_ActiveRotationMatrix<DataType>& operator/=(const ASD_ActiveRotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)/= o.data(i); return *this; }
    ASD_ActiveRotationMatrix<DataType>& operator*=(const ASD_ActiveRotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)*= o.data(i); return *this; }
    ASD_ActiveRotationMatrix<DataType> operator/(const ASD_ActiveRotationMatrix<DataType>& o) const { ASD_ActiveRotationMatrix<DataType> out(*this); return out /= o; }
    ASD_ActiveRotationMatrix<DataType> operator*(const ASD_ActiveRotationMatrix<DataType>& o) const { ASD_ActiveRotationMatrix<DataType> out(*this); return out *= o; }
    ASD_ActiveRotationMatrix<DataType>& operator=(const ASD_ActiveRotationMatrix<DataType>& o) { std::copy_n(o.data(), size(), data());  return *this; }

    // size of the file
    int size() const { return size_; }
    // zero out
    void zero() { fill(DataType(0.0)); }
    void fill(const DataType& a) { std::fill_n(data(), size_, a); }
    // returns dot product
    DataType dot_product(const ASD_ActiveRotationMatrix<DataType>& o) const { return blas::dot_product(data(), size_, o.data()); }
    DataType dot_product(const std::shared_ptr<const ASD_ActiveRotationMatrix<DataType>> o) const { return dot_product(*o); }
    // scale function
    void scale(const DataType& a) { std::for_each(data(), data()+size_, [&a](DataType& p) { p *= a; }); }
    // returns norm of the vector
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    // returns a root mean square
    double rms() const { return norm() / std::sqrt(static_cast<double>(size())); }
    // daxpy added to self
    void ax_plus_y(const DataType& a, const ASD_ActiveRotationMatrix& o) { std::transform(o.data(), o.data()+size_, data(), data(), [&a](DataType p, DataType q) { return p*a+q; }); }
    void ax_plus_y(const DataType& a, const std::shared_ptr<const ASD_ActiveRotationMatrix> o) { ax_plus_y(a, *o); }

    // orthogonalize to the liset of ASD_ActiveRotationMatrix's
    double orthog(std::list<std::shared_ptr<const ASD_ActiveRotationMatrix<DataType>>> c) {
      for (auto iter = c.begin(); iter != c.end(); ++iter)
        this->ax_plus_y(- this->dot_product(**iter), **iter);
      return normalize();
    }

    double normalize() {
      const double scal = 1.0/this->norm();
      scale(scal);
      return 1.0/scal;
    }

    void synchronize() {
#ifdef HAVE_MPI_H
      mpi__->broadcast(data(), size(), 0);
#endif
    }

    // return data_
    DataType* data() { return data_.get(); }
    const DataType* data() const { return data_.get(); }
    DataType& data(const size_t i) { return data_[i]; }
    const DataType& data(const size_t i) const { return data_[i]; }
    // return data_
    DataType* begin() { return data(); }
    // return data_
    DataType* end() { return data()+size_; }

    // ROTATION MATRIX
    //   C   A   V        A(A) A(B)   ORDERED AS:
    // C .  ca  cv    A(A)  .  [aa]   ca -> va -> vc -> aa
    // A . [aa] av    A(B)  .    .
    // V .   .   .

    // active-active block. A index runs first
    DataType* ptr_aa() { return data(); }
    DataType& ele_aa(const int ib, const int ia) { return data_[ib + ia*nactB_]; }

    // const references and pointers
    const DataType& ele_aa(const int ib, const int ia) const { return data_[ib + ia*nactB_]; }
    const DataType* ptr_aa() const { return data(); }

    // unpack to Matrix
    template<class MatType>
    std::shared_ptr<MatType> unpack(const DataType a = 0.0) const {
    //const int nocc = nclosed_ + nact_;
      const int nbasis = nclosed_ + nact_ + nvirt_;
      auto out = std::make_shared<MatType>(nbasis, nbasis);
      std::fill_n(out->data(), out->size(), a);

      //Active-Active 
      for (int i = 0; i != nactA_; ++i) {
        for (int j = 0; j != nactB_; ++j) {
          out->element(j+nclosed_+nactA_,i+nclosed_) = ele_aa(j,i);
        }
      }
      //Anti-symmetric
      for (int i = 0; i != nbasis; ++i) {
        for (int j = 0; j <= i; ++j) {
          out->element(j, i) = - detail::conj(out->element(i, j));
        }
      }

      return out;
    }

    // print matrix
    void print(const std::string in = "") const {
      std::cout << "++++ " + in + " ++++" << std::endl;
      if (nactA_ && nactB_) {
        std::cout << " printing active(B)-active(A) block" << std::endl;
        for (int i = 0; i != nactA_; ++i) {
          for (int j = 0; j != nactB_; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa(j,i);
          }
          std::cout << std::endl;
        }
      }

    }
};

using ASD_ActiveRotFile = ASD_ActiveRotationMatrix<double>;

}

#endif
