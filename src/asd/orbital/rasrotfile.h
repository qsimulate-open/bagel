//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbital/rasrotfile.h
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

#ifndef __BAGEL_ASD_RASROTFILE_H
#define __BAGEL_ASD_RASROTFILE_H

#include <src/util/math/matrix.h>

namespace bagel {

template<typename DataType>
class ASD_RAS_RotationMatrix {
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

    std::array<int,3> rasA_;
    std::array<int,3> rasB_;

  public:
    ASD_RAS_RotationMatrix(const int iclos, const int iact, const int ivirt, const int iactA, const int iactB, std::array<int,3> rasA, std::array<int,3> rasB)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), nactA_(iactA), nactB_(iactB), size_(iclos*iact+iclos*ivirt+iact*ivirt+iactA*iactB 
       + rasA[0]*rasA[1] + rasA[0]*rasA[2] + rasA[1]*rasA[2]
       + rasB[0]*rasB[1] + rasB[0]*rasB[2] + rasB[1]*rasB[2]), data_(new DataType[size_]), rasA_(rasA), rasB_(rasB) {
      zero();
    }
    ASD_RAS_RotationMatrix(const ASD_RAS_RotationMatrix& o) : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), nactA_(o.nactA_), nactB_(o.nactB_),  size_(o.size_), data_(new DataType[o.size_]), rasA_(o.rasA_), rasB_(o.rasB_) {
      *this = o;
    }
    ASD_RAS_RotationMatrix(std::shared_ptr<const ASD_RAS_RotationMatrix> o)
      : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), nactA_(o->nactA_), nactB_(o->nactB_), size_(o->size_), data_(new DataType[o->size_]), rasA_(o->rasA_), rasB_(o->rasB_) {
      *this = *o;
    }
    //Matrix to ASDRotationMatrix conversion
    ASD_RAS_RotationMatrix(std::shared_ptr<const Matrix_base<DataType>> o, const int iclos, const int iact, const int ivirt, const int iactA, const int iactB, std::array<int,3> rasA, std::array<int,3> rasB)
      : nclosed_(iclos), nact_(iact), nvirt_(ivirt), nactA_(iactA), nactB_(iactB), size_(iclos*iact+iclos*ivirt+iact*ivirt+iactA*iactB
        + rasA[0]*rasA[1] + rasA[0]*rasA[2] + rasA[1]*rasA[2] 
        + rasB[0]*rasB[1] + rasB[0]*rasB[2] + rasB[1]*rasB[2]), data_(new DataType[size_]), rasA_(rasA), rasB_(rasB) {

      const int nocc = nclosed_ + nact_;
      assert(iact == iactA+iactB);
      for (int i = 0; i != nact_; ++i) {
        for (int j = 0; j != nvirt_;   ++j) {
          ele_va(j, i) = o->element(j+nocc, i+nclosed_);
        }
        for (int j = 0; j != nclosed_; ++j) {
          ele_ca(j, i) = o->element(i+nclosed_, j);
        }
      }
      for (int i = 0; i != nclosed_; ++i) {
        for (int j = 0; j != nvirt_;   ++j) {
          ele_vc(j, i) = o->element(j+nocc, i);
        }
      }
      for (int i = 0; i != nactA_; ++i) {
        for (int j = 0; j != nactB_; ++j) {
          ele_aa(j, i) = o->element(j+nclosed_+nactA_,i+nclosed_);
        }
      }
      //TODO: template aaXX 12, 13, 23 with <0> for A <1> for B
      //A
      for (int i = 0; i != rasA_[0]; ++i) 
        for (int j = 0; j != rasA_[1]; ++ j) {
          ele_aa21A(j, i) = o->element(j+nclosed_+rasA_[0], i+nclosed_);
        }
      for (int i = 0; i != rasA_[0]; ++i) 
        for (int j = 0; j != rasA_[2]; ++ j) {
          ele_aa31A(j, i) = o->element(j+nclosed_+rasA_[0]+rasA_[1], i+nclosed_);
        }
      for (int i = 0; i != rasA_[1]; ++i) 
        for (int j = 0; j != rasA_[2]; ++ j) {
          ele_aa32A(j, i) = o->element(j+nclosed_+rasA_[0]+rasA_[1], i+nclosed_+rasA_[0]);
        }
      //B
      for (int i = 0; i != rasB_[0]; ++i) 
        for (int j = 0; j != rasB_[1]; ++ j) {
          ele_aa21B(j, i) = o->element(j+nclosed_+nactA_+rasB_[0], i+nclosed_+nactA_);
        }
      for (int i = 0; i != rasB_[0]; ++i) 
        for (int j = 0; j != rasB_[2]; ++ j) {
          ele_aa31B(j, i) = o->element(j+nclosed_+nactA_+rasB_[0]+rasB_[1], i+nclosed_+nactA_);
        }
      for (int i = 0; i != rasB_[1]; ++i) 
        for (int j = 0; j != rasB_[2]; ++ j) {
          ele_aa32B(j, i) = o->element(j+nclosed_+nactA_+rasB_[0]+rasB_[1], i+nclosed_+nactA_+rasB_[0]);
        }
    }

    std::shared_ptr<ASD_RAS_RotationMatrix<DataType>> clone() const {
      return std::make_shared<ASD_RAS_RotationMatrix<DataType>>(nclosed_, nact_, nvirt_, nactA_, nactB_, rasA_, rasB_);
    }
    std::shared_ptr<ASD_RAS_RotationMatrix<DataType>> copy() const {
      return std::make_shared<ASD_RAS_RotationMatrix<DataType>>(*this);
    }

    // overloaded operators
    ASD_RAS_RotationMatrix<DataType> operator+(const ASD_RAS_RotationMatrix<DataType>& o) const { ASD_RAS_RotationMatrix<DataType> out(*this); out.ax_plus_y(1.0, o); return out; }
    ASD_RAS_RotationMatrix<DataType> operator-(const ASD_RAS_RotationMatrix<DataType>& o) const { ASD_RAS_RotationMatrix<DataType> out(*this); out.ax_plus_y(-1.0, o); return out; }
    ASD_RAS_RotationMatrix<DataType>& operator+=(const ASD_RAS_RotationMatrix<DataType>& o) { ax_plus_y(1.0, o); return *this; }
    ASD_RAS_RotationMatrix<DataType>& operator-=(const ASD_RAS_RotationMatrix<DataType>& o) { ax_plus_y(-1.0, o); return *this; }
    ASD_RAS_RotationMatrix<DataType>& operator*=(const DataType a) { scale(a); return *this; }
    ASD_RAS_RotationMatrix<DataType>& operator/=(const ASD_RAS_RotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)/= o.data(i); return *this; }
    ASD_RAS_RotationMatrix<DataType>& operator*=(const ASD_RAS_RotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)*= o.data(i); return *this; }
    ASD_RAS_RotationMatrix<DataType> operator/(const ASD_RAS_RotationMatrix<DataType>& o) const { ASD_RAS_RotationMatrix<DataType> out(*this); return out /= o; }
    ASD_RAS_RotationMatrix<DataType> operator*(const ASD_RAS_RotationMatrix<DataType>& o) const { ASD_RAS_RotationMatrix<DataType> out(*this); return out *= o; }
    ASD_RAS_RotationMatrix<DataType>& operator=(const ASD_RAS_RotationMatrix<DataType>& o) { std::copy_n(o.data(), size(), data());  return *this; }

    // size of the file
    int size() const { return size_; }
    // zero out
    void zero() { fill(DataType(0.0)); }
    void fill(const DataType& a) { std::fill_n(data(), size_, a); }
    // returns dot product
    DataType dot_product(const ASD_RAS_RotationMatrix<DataType>& o) const { return blas::dot_product(data(), size_, o.data()); }
    DataType dot_product(const std::shared_ptr<const ASD_RAS_RotationMatrix<DataType>> o) const { return dot_product(*o); }
    // scale function
    void scale(const DataType& a) { std::for_each(data(), data()+size_, [&a](DataType& p) { p *= a; }); }
    // returns norm of the vector
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    // returns a root mean square
    double rms() const { return norm() / std::sqrt(static_cast<double>(size())); }
    // daxpy added to self
    void ax_plus_y(const DataType& a, const ASD_RAS_RotationMatrix& o) { std::transform(o.data(), o.data()+size_, data(), data(), [&a](DataType p, DataType q) { return p*a+q; }); }
    void ax_plus_y(const DataType& a, const std::shared_ptr<const ASD_RAS_RotationMatrix> o) { ax_plus_y(a, *o); }

    // orthogonalize to the liset of ASD_RAS_RotationMatrix's
    double orthog(std::list<std::shared_ptr<const ASD_RAS_RotationMatrix<DataType>>> c) {
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

    // closed-active block. closed runs first
    DataType* ptr_ca() { return data(); }
    DataType& ele_ca(const int ic, const int ia) { return data_[ic + ia*nclosed_]; }
    // active-virtual block. virtual runs first
    DataType* ptr_va() { return data() + nclosed_*nact_; }
    DataType& ele_va(const int iv, const int ia) { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    // closed-virtual block. virtual runs first
    DataType* ptr_vc() { return data() + (nclosed_+nvirt_)*nact_; }
    DataType& ele_vc(const int iv, const int ic) { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    // active-active block. A index runs first
    DataType* ptr_aa() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_; }
    DataType& ele_aa(const int ib, const int ia) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ib + ia*nactB_]; }
    // RAS A
    DataType* ptr_aa21A() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_; }
    DataType& ele_aa21A(const int i2, const int i1) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ 
                                                                  + i2 + i1*rasA_[1]]; }

    DataType* ptr_aa31A() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1]; }
    DataType& ele_aa31A(const int i3, const int i1) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ 
                                                                  + rasA_[0]*rasA_[1] + i3 + i1*rasA_[2]]; }

    DataType* ptr_aa32A() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2]; }
    DataType& ele_aa32A(const int i3, const int i2) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ 
                                                                  + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + i3 + i2*rasA_[2]]; }
    // RASB
    DataType* ptr_aa21B() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]; }
    DataType& ele_aa21B(const int i2, const int i1) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ 
                                                                  + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] 
                                                                  + i2 + i1*rasB_[1]]; }

    DataType* ptr_aa31B() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]
                                                                                                      + rasB_[0]*rasB_[1]; }
    DataType& ele_aa31B(const int i3, const int i1) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ 
                                                                  + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] 
                                                                  + rasB_[0]*rasB_[1] + i3 + i1*rasB_[2]]; }

    DataType* ptr_aa32B() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]
                                                                                                      + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2]; }
    DataType& ele_aa32B(const int i3, const int i2) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ 
                                                                  + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] 
                                                                  + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2] + i3 + i2*rasB_[2]]; }


    // const references and pointers
    const DataType& ele_ca(const int ic, const int ia) const { return data_[ic + ia*nclosed_]; }
    const DataType& ele_va(const int iv, const int ia) const { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    const DataType& ele_vc(const int iv, const int ic) const { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    const DataType& ele_aa(const int ib, const int ia) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ib + ia*nactB_]; }

    const DataType& ele_aa21A(const int i2, const int i1) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + i2 + i1*rasA_[1]]; }
    const DataType& ele_aa31A(const int i3, const int i1) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + i3 + i1*rasA_[2]]; }
    const DataType& ele_aa32A(const int i3, const int i2) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + i3 + i2*rasA_[2]]; }
    const DataType& ele_aa21B(const int i2, const int i1) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + i2 + i1*rasB_[1]]; }
    const DataType& ele_aa31B(const int i3, const int i1) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1] + i3 + i1*rasB_[2]]; }
    const DataType& ele_aa32B(const int i3, const int i2) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2] + i3 + i2*rasB_[2]]; }

    const DataType* ptr_ca() const { return data(); }
    const DataType* ptr_va() const { return data() + nclosed_*nact_; }
    const DataType* ptr_vc() const { return data() + (nclosed_+nvirt_)*nact_; }
    const DataType* ptr_aa() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_; }
    const DataType* ptr_aa21A() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_; }
    const DataType* ptr_aa31A() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1]; }
    const DataType* ptr_aa32A() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2]; }
    const DataType* ptr_aa21B() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]; }
    const DataType* ptr_aa31B() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]
                                                                                                                  + rasB_[0]*rasB_[1]; }
    const DataType* ptr_aa32B() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]
                                                                                                                  + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2]; }

    // unpack to Matrix
    template<class MatType>
    std::shared_ptr<MatType> unpack(const DataType a = 0.0) const {
      const int nocc = nclosed_ + nact_;
      const int nbasis = nclosed_ + nact_ + nvirt_;
      auto out = std::make_shared<MatType>(nbasis, nbasis);
      std::fill_n(out->data(), out->size(), a);

      for (int i = 0; i != nact_; ++i) {
        //Virtual-Active
        for (int j = 0; j != nvirt_;   ++j) {
          out->element(j+nocc, i+nclosed_) = ele_va(j, i);
        }
        //Active-Closed
        for (int j = 0; j != nclosed_; ++j) {
          out->element(i+nclosed_, j) = ele_ca(j, i);
        }
      }
      //Virtual-Closed
      for (int i = 0; i != nclosed_; ++i) {
        for (int j = 0; j != nvirt_;   ++j) {
          out->element(j+nocc, i) = ele_vc(j, i);
        }
      }
      //Active-Active 
      for (int i = 0; i != nactA_; ++i) {
        for (int j = 0; j != nactB_; ++j) {
          out->element(j+nclosed_+nactA_,i+nclosed_) = ele_aa(j,i);
        }
      }
      //RAS A
      for (int i = 0; i != rasA_[0]; ++i)
        for (int j = 0; j != rasA_[1]; ++j) {
          out->element(j+nclosed_+rasA_[0],i+nclosed_) = ele_aa21A(j,i);
        }
      for (int i = 0; i != rasA_[0]; ++i) 
        for (int j = 0; j != rasA_[2]; ++j) {
          out->element(j+nclosed_+rasA_[0]+rasA_[1],i+nclosed_) = ele_aa31A(j,i);
        }
      for (int i = 0; i != rasA_[1]; ++i) 
        for (int j = 0; j != rasA_[2]; ++j) {
          out->element(j+nclosed_+rasA_[0]+rasA_[1],i+nclosed_+rasA_[0]) = ele_aa32A(j,i);
        }
      //RAS B
      for (int i = 0; i != rasB_[0]; ++i) 
        for (int j = 0; j != rasB_[1]; ++j) {
          out->element(j+nclosed_+nactA_+rasB_[0],i+nclosed_+nactA_) = ele_aa21B(j,i);
        }
      for (int i = 0; i != rasB_[0]; ++i) 
        for (int j = 0; j != rasB_[2]; ++j) {
          out->element(j+nclosed_+nactA_+rasB_[0]+rasB_[1],i+nclosed_+nactA_) = ele_aa31B(j,i);
        }
      for (int i = 0; i != rasB_[1]; ++i) 
        for (int j = 0; j != rasB_[2]; ++j) {
          out->element(j+nclosed_+nactA_+rasB_[0]+rasB_[1],i+nclosed_+nactA_+rasB_[0]) = ele_aa32B(j,i);
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
};

using ASD_RAS_RotFile = ASD_RAS_RotationMatrix<double>;

}

#endif
