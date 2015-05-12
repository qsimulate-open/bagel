//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbital/asd_rotfile.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim: <inkoo.kim@northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifndef __BAGEL_ASD_ROTFILE_H
#define __BAGEL_ASD_ROTFILE_H

#include <src/util/math/matrix.h>

namespace bagel {

template<typename DataType>
class ASD_RotationMatrix {
  public:
    using data_type = DataType;

  protected:
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    std::array<int,3> rasA_;
    std::array<int,3> rasB_;
    const int nactA_;
    const int nactB_;
    const bool inter_; //inter subspace rotations
    const bool intra_; //intra active space rotations
    const int size_;
    std::unique_ptr<DataType[]> data_;


  public:
    ASD_RotationMatrix(const int iclos, const int iact, const int ivirt, std::array<int,3> rasA, std::array<int,3> rasB, const bool inter, const bool intra)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), rasA_(rasA), rasB_(rasB), nactA_(rasA[0]+rasA[1]+rasA[2]), nactB_(rasB[0]+rasB[1]+rasB[2]), inter_(inter), intra_(intra),
       size_((inter_ ? iclos*iact + iclos*ivirt + iact*ivirt : 0) + (intra_ ? nactA_*nactB_ + rasA[0]*rasA[1] + rasA[0]*rasA[2] + rasA[1]*rasA[2] + rasB[0]*rasB[1] + rasB[0]*rasB[2] + rasB[1]*rasB[2] : 0)), data_(new DataType[size_]) {
      assert(!(!inter_ && !intra_));
      assert(nact_ == rasA[0]+rasA[1]+rasA[2]+rasB[0]+rasB[1]+rasB[2]);
      assert(nact_ == nactA_ + nactB_);
      zero();
    }
    ASD_RotationMatrix(const ASD_RotationMatrix& o)
      : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), rasA_(o.rasA_), rasB_(o.rasB_), nactA_(o.nactA_), nactB_(o.nactB_), inter_(o.inter_), intra_(o.intra_), size_(o.size_), data_(new DataType[o.size_]) {
      *this = o;
    }
    ASD_RotationMatrix(std::shared_ptr<const ASD_RotationMatrix> o)
      : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), rasA_(o->rasA_), rasB_(o->rasB_), nactA_(o->nactA_), nactB_(o->nactB_), inter_(o->inter_), intra_(o->intra_), size_(o->size_), data_(new DataType[o->size_]) {
      *this = *o;
    }
    //Matrix to ASDRotationMatrix conversion
    ASD_RotationMatrix(std::shared_ptr<const Matrix_base<DataType>> o, const int iclos, const int iact, const int ivirt, std::array<int,3> rasA, std::array<int,3> rasB, const bool inter, const bool intra)
      : nclosed_(iclos), nact_(iact), nvirt_(ivirt), rasA_(rasA), rasB_(rasB), nactA_(rasA[0]+rasA[1]+rasA[2]), nactB_(rasB[0]+rasB[1]+rasB[2]), inter_(inter), intra_(intra),
        size_((inter_ ? iclos*iact + iclos*ivirt + iact*ivirt : 0) + (intra_ ? nactA_*nactB_ + rasA[0]*rasA[1] + rasA[0]*rasA[2] + rasA[1]*rasA[2] + rasB[0]*rasB[1] + rasB[0]*rasB[2] + rasB[1]*rasB[2] : 0)), data_(new DataType[size_]) {
      assert(!(!inter_ && !intra_));
      assert(nact_ == rasA[0]+rasA[1]+rasA[2]+rasB[0]+rasB[1]+rasB[2]);
      assert(nact_ == nactA_ + nactB_);
      const int nocc = nclosed_ + nact_;
      if (inter_) {
        for (int i = 0; i != nact_; ++i) {
          for (int j = 0; j != nvirt_;   ++j)
            ele_va(j, i) = o->element(j+nocc, i+nclosed_);
          for (int j = 0; j != nclosed_; ++j)
            ele_ca(j, i) = o->element(i+nclosed_, j);
        }
        for (int i = 0; i != nclosed_; ++i)
          for (int j = 0; j != nvirt_;   ++j)
            ele_vc(j, i) = o->element(j+nocc, i);
      }

      if (intra_) {
        for (int i = 0; i != nactA_; ++i)
          for (int j = 0; j != nactB_; ++j)
            ele_aa(j, i) = o->element(j+nclosed_+nactA_,i+nclosed_);
        //A
        for (int i = 0; i != rasA_[0]; ++i)
          for (int j = 0; j != rasA_[1]; ++ j)
            ele_aa21A(j, i) = o->element(j+nclosed_+rasA_[0], i+nclosed_);
        for (int i = 0; i != rasA_[0]; ++i)
          for (int j = 0; j != rasA_[2]; ++ j)
            ele_aa31A(j, i) = o->element(j+nclosed_+rasA_[0]+rasA_[1], i+nclosed_);
        for (int i = 0; i != rasA_[1]; ++i)
          for (int j = 0; j != rasA_[2]; ++ j)
            ele_aa32A(j, i) = o->element(j+nclosed_+rasA_[0]+rasA_[1], i+nclosed_+rasA_[0]);
        //B
        for (int i = 0; i != rasB_[0]; ++i)
          for (int j = 0; j != rasB_[1]; ++ j)
            ele_aa21B(j, i) = o->element(j+nclosed_+nactA_+rasB_[0], i+nclosed_+nactA_);
        for (int i = 0; i != rasB_[0]; ++i)
          for (int j = 0; j != rasB_[2]; ++ j)
            ele_aa31B(j, i) = o->element(j+nclosed_+nactA_+rasB_[0]+rasB_[1], i+nclosed_+nactA_);
        for (int i = 0; i != rasB_[1]; ++i)
          for (int j = 0; j != rasB_[2]; ++ j)
            ele_aa32B(j, i) = o->element(j+nclosed_+nactA_+rasB_[0]+rasB_[1], i+nclosed_+nactA_+rasB_[0]);
      }
    }

    std::shared_ptr<ASD_RotationMatrix<DataType>> clone() const {
      return std::make_shared<ASD_RotationMatrix<DataType>>(nclosed_, nact_, nvirt_, rasA_, rasB_, inter_, intra_);
    }
    std::shared_ptr<ASD_RotationMatrix<DataType>> copy() const {
      return std::make_shared<ASD_RotationMatrix<DataType>>(*this);
    }

    // overloaded operators
    ASD_RotationMatrix<DataType> operator+(const ASD_RotationMatrix<DataType>& o) const { ASD_RotationMatrix<DataType> out(*this); out.ax_plus_y(1.0, o); return out; }
    ASD_RotationMatrix<DataType> operator-(const ASD_RotationMatrix<DataType>& o) const { ASD_RotationMatrix<DataType> out(*this); out.ax_plus_y(-1.0, o); return out; }
    ASD_RotationMatrix<DataType>& operator+=(const ASD_RotationMatrix<DataType>& o) { ax_plus_y(1.0, o); return *this; }
    ASD_RotationMatrix<DataType>& operator-=(const ASD_RotationMatrix<DataType>& o) { ax_plus_y(-1.0, o); return *this; }
    ASD_RotationMatrix<DataType>& operator*=(const DataType a) { scale(a); return *this; }
    ASD_RotationMatrix<DataType>& operator/=(const ASD_RotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)/= o.data(i); return *this; }
    ASD_RotationMatrix<DataType>& operator*=(const ASD_RotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)*= o.data(i); return *this; }
    ASD_RotationMatrix<DataType> operator/(const ASD_RotationMatrix<DataType>& o) const { ASD_RotationMatrix<DataType> out(*this); return out /= o; }
    ASD_RotationMatrix<DataType> operator*(const ASD_RotationMatrix<DataType>& o) const { ASD_RotationMatrix<DataType> out(*this); return out *= o; }
    ASD_RotationMatrix<DataType>& operator=(const ASD_RotationMatrix<DataType>& o) { std::copy_n(o.data(), size(), data());  return *this; }

    // size of the file
    int size() const { return size_; }
    // zero out
    void zero() { fill(DataType(0.0)); }
    void fill(const DataType& a) { std::fill_n(data(), size_, a); }
    // returns dot product
    DataType dot_product(const ASD_RotationMatrix<DataType>& o) const { return blas::dot_product(data(), size_, o.data()); }
    DataType dot_product(const std::shared_ptr<const ASD_RotationMatrix<DataType>> o) const { return dot_product(*o); }
    // scale function
    void scale(const DataType& a) { std::for_each(data(), data()+size_, [&a](DataType& p) { p *= a; }); }
    // returns norm of the vector
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    // returns a root mean square
    double rms() const { return norm() / std::sqrt(static_cast<double>(size())); }
    // daxpy added to self
    void ax_plus_y(const DataType& a, const ASD_RotationMatrix& o) { std::transform(o.data(), o.data()+size_, data(), data(), [&a](DataType p, DataType q) { return p*a+q; }); }
    void ax_plus_y(const DataType& a, const std::shared_ptr<const ASD_RotationMatrix> o) { ax_plus_y(a, *o); }

    // orthogonalize to the liset of ASD_RotationMatrix's
    double orthog(std::list<std::shared_ptr<const ASD_RotationMatrix<DataType>>> c) {
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

    bool is_inter() const { return inter_; }
    bool is_intra() const { return intra_; }

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
    DataType* ptr_ca() { assert(inter_); return data(); }
    DataType& ele_ca(const int ic, const int ia) { assert(inter_); return data_[ic + ia*nclosed_]; }
    // active-virtual block. virtual runs first
    DataType* ptr_va() { assert(inter_); return data() + nclosed_*nact_; }
    DataType& ele_va(const int iv, const int ia) { assert(inter_); return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    // closed-virtual block. virtual runs first
    DataType* ptr_vc() { assert(inter_); return data() + (nclosed_+nvirt_)*nact_; }
    DataType& ele_vc(const int iv, const int ic) { assert(inter_); return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }

    // active-active block. B index runs first
    DataType* ptr_aa() { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0); }
    DataType& ele_aa(const int ib, const int ia) { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + ib + ia*nactB_]; }
    // RAS A
    DataType* ptr_aa21A() { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_; }
    DataType* ptr_aa31A() { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1]; }
    DataType* ptr_aa32A() { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2]; }

    DataType& ele_aa21A(const int i2, const int i1) { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + i2 + i1*rasA_[1]]; }
    DataType& ele_aa31A(const int i3, const int i1) { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + i3 + i1*rasA_[2]]; }
    DataType& ele_aa32A(const int i3, const int i2) { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + i3 + i2*rasA_[2]]; }
    // RASB
    DataType* ptr_aa21B() { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]; }
    DataType* ptr_aa31B() { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1]; }
    DataType* ptr_aa32B() { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2]; }

    DataType& ele_aa21B(const int i2, const int i1) { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + i2 + i1*rasB_[1]]; }
    DataType& ele_aa31B(const int i3, const int i1) { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1] + i3 + i1*rasB_[2]]; }
    DataType& ele_aa32B(const int i3, const int i2) { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2] + i3 + i2*rasB_[2]]; }

    // const references and pointer_s
    const DataType& ele_ca(const int ic, const int ia) const { assert(inter_); return data_[ic + ia*nclosed_]; }
    const DataType& ele_va(const int iv, const int ia) const { assert(inter_); return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    const DataType& ele_vc(const int iv, const int ic) const { assert(inter_); return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }

    const DataType& ele_aa   (const int ib, const int ia) const { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + ib + ia*nactB_]; }
    const DataType& ele_aa21A(const int i2, const int i1) const { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + i2 + i1*rasA_[1]]; }
    const DataType& ele_aa31A(const int i3, const int i1) const { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + i3 + i1*rasA_[2]]; }
    const DataType& ele_aa32A(const int i3, const int i2) const { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + i3 + i2*rasA_[2]]; }
    const DataType& ele_aa21B(const int i2, const int i1) const { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + i2 + i1*rasB_[1]]; }
    const DataType& ele_aa31B(const int i3, const int i1) const { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1] + i3 + i1*rasB_[2]]; }
    const DataType& ele_aa32B(const int i3, const int i2) const { assert(intra_); return data_[(inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2] + i3 + i2*rasB_[2]]; }

    const DataType* ptr_ca() const { assert(inter_); return data(); }
    const DataType* ptr_va() const { assert(inter_); return data() + nclosed_*nact_; }
    const DataType* ptr_vc() const { assert(inter_); return data() + (nclosed_+nvirt_)*nact_; }

    const DataType* ptr_aa()    const { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0); }
    const DataType* ptr_aa21A() const { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_; }
    const DataType* ptr_aa31A() const { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1]; }
    const DataType* ptr_aa32A() const { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2]; }
    const DataType* ptr_aa21B() const { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]; }
    const DataType* ptr_aa31B() const { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1]; }
    const DataType* ptr_aa32B() const { assert(intra_); return data() + (inter_ ? (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ : 0) + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2]; }

    // unpack to Matrix
    template<class MatType>
    std::shared_ptr<MatType> unpack(const DataType a = 0.0) const {
      const int nocc = nclosed_ + nact_;
      const int nbasis = nclosed_ + nact_ + nvirt_;
      auto out = std::make_shared<MatType>(nbasis, nbasis);
      std::fill_n(out->data(), out->size(), a);

      if (inter_) {
        for (int i = 0; i != nact_; ++i) {
          for (int j = 0; j != nvirt_;   ++j)
            out->element(j+nocc, i+nclosed_) = ele_va(j, i); //Virtual-Active
          for (int j = 0; j != nclosed_; ++j)
            out->element(i+nclosed_, j) = ele_ca(j, i); //Active-Closed
        }
        for (int i = 0; i != nclosed_; ++i)
          for (int j = 0; j != nvirt_;   ++j)
            out->element(j+nocc, i) = ele_vc(j, i); //Virtual-Closed
      }

      if (intra_) {
        for (int i = 0; i != nactA_; ++i)
          for (int j = 0; j != nactB_; ++j)
            out->element(j+nclosed_+nactA_,i+nclosed_) = ele_aa(j,i); //Active-Active
        //RAS A
        for (int i = 0; i != rasA_[0]; ++i)
          for (int j = 0; j != rasA_[1]; ++j)
            out->element(j+nclosed_+rasA_[0],i+nclosed_) = ele_aa21A(j,i);
        for (int i = 0; i != rasA_[0]; ++i)
          for (int j = 0; j != rasA_[2]; ++j)
            out->element(j+nclosed_+rasA_[0]+rasA_[1],i+nclosed_) = ele_aa31A(j,i);
        for (int i = 0; i != rasA_[1]; ++i)
          for (int j = 0; j != rasA_[2]; ++j)
            out->element(j+nclosed_+rasA_[0]+rasA_[1],i+nclosed_+rasA_[0]) = ele_aa32A(j,i);
        //RAS B
        for (int i = 0; i != rasB_[0]; ++i)
          for (int j = 0; j != rasB_[1]; ++j)
            out->element(j+nclosed_+nactA_+rasB_[0],i+nclosed_+nactA_) = ele_aa21B(j,i);
        for (int i = 0; i != rasB_[0]; ++i)
          for (int j = 0; j != rasB_[2]; ++j)
            out->element(j+nclosed_+nactA_+rasB_[0]+rasB_[1],i+nclosed_+nactA_) = ele_aa31B(j,i);
        for (int i = 0; i != rasB_[1]; ++i)
          for (int j = 0; j != rasB_[2]; ++j)
            out->element(j+nclosed_+nactA_+rasB_[0]+rasB_[1],i+nclosed_+nactA_+rasB_[0]) = ele_aa32B(j,i);
      }

      //Anti-symmetric
      for (int i = 0; i != nbasis; ++i)
        for (int j = 0; j <= i; ++j)
          out->element(j, i) = - detail::conj(out->element(i, j));

      return out;
    }

    // print matrix
    void print(const std::string in = "") const {
      std::cout << "++++ " + in + " ++++" << std::endl;
      if (inter_) {
        if (nact_ && nclosed_) {
          std::cout << " printing closed-active block" << std::endl;
          for (int i = 0; i != nact_; ++i) {
            for (int j = 0; j != nclosed_; ++j) {
              std::cout << std::setw(10) << std::setprecision(6) << ele_ca(j,i);
            }
            std::cout << std::endl;
          }
        }
        if (nact_ && nvirt_) {
          std::cout << " printing virtual-active block" << std::endl;
          for (int i = 0; i != nact_; ++i) {
            for (int j = 0; j != nvirt_; ++j) {
              std::cout << std::setw(10) << std::setprecision(6) << ele_va(j,i);
            }
            std::cout << std::endl;
          }
        }
        if (nclosed_ && nvirt_) {
          std::cout << " printing virtual-closed block" << std::endl;
          for (int i = 0; i != nclosed_; ++i) {
            for (int j = 0; j != nvirt_; ++j) {
              std::cout << std::setw(10) << std::setprecision(6) << ele_vc(j,i);
            }
            std::cout << std::endl;
          }
        }
      }

      if (intra_) {
        if (nactA_ && nactB_) {
          std::cout << " printing active(B)-active(A) block" << std::endl;
          for (int i = 0; i != nactA_; ++i) {
            for (int j = 0; j != nactB_; ++j) {
              std::cout << std::setw(10) << std::setprecision(6) << ele_aa(j,i);
            }
            std::cout << std::endl;
          }
        }
        if (nactA_ && nactB_) {
          std::cout << " printing RAS(A) blocks" << std::endl;
          std::cout << "          RAS1-RAS2 blocks" << std::endl;
          for (int i = 0; i != rasA_[0]; ++i)
            for (int j = 0; j != rasA_[1]; ++j) {
              std::cout << std::setw(10) << std::setprecision(6) << ele_aa21A(j,i);
            }
          std::cout << std::endl;
          std::cout << "          RAS1-RAS3 blocks" << std::endl;
          for (int i = 0; i != rasA_[0]; ++i)
            for (int j = 0; j != rasA_[2]; ++j) {
              std::cout << std::setw(10) << std::setprecision(6) << ele_aa31A(j,i);
            }
          std::cout << std::endl;
          std::cout << "          RAS2-RAS3 blocks" << std::endl;
          for (int i = 0; i != rasA_[1]; ++i)
            for (int j = 0; j != rasA_[2]; ++j) {
              std::cout << std::setw(10) << std::setprecision(6) << ele_aa32A(j,i);
            }
          std::cout << std::endl;
          std::cout << " printing RAS(B) blocks" << std::endl;
          std::cout << "          RAS1-RAS2 blocks" << std::endl;
          for (int i = 0; i != rasB_[0]; ++i)
            for (int j = 0; j != rasB_[1]; ++j) {
              std::cout << std::setw(10) << std::setprecision(6) << ele_aa21B(j,i);
            }
          std::cout << std::endl;
          std::cout << "          RAS1-RAS3 blocks" << std::endl;
          for (int i = 0; i != rasB_[0]; ++i)
            for (int j = 0; j != rasB_[2]; ++j) {
              std::cout << std::setw(10) << std::setprecision(6) << ele_aa31B(j,i);
            }
          std::cout << std::endl;
          std::cout << "          RAS2-RAS3 blocks" << std::endl;
          for (int i = 0; i != rasB_[1]; ++i)
            for (int j = 0; j != rasB_[2]; ++j) {
              std::cout << std::setw(10) << std::setprecision(6) << ele_aa32B(j,i);
            }
          std::cout << std::endl;
        }
      }

    }
};

using ASD_RotFile = ASD_RotationMatrix<double>;

}

#endif
