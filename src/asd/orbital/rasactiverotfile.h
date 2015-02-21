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

#ifndef __BAGEL_ASD_RAS_ACTIVEROTFILE_H
#define __BAGEL_ASD_RAS_ACTIVEROTFILE_H

#include <src/util/math/matrix.h>

namespace bagel {

template<typename DataType>
class ASD_RAS_ActiveRotationMatrix {
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
    ASD_RAS_ActiveRotationMatrix(const int iclos, const int iact, const int ivirt, const int iactA, const int iactB, std::array<int,3> rasA, std::array<int,3> rasB)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), nactA_(iactA), nactB_(iactB), size_(iactA*iactB 
       + rasA[0]*rasA[1] + rasA[0]*rasA[2] + rasA[1]*rasA[2]
       + rasB[0]*rasB[1] + rasB[0]*rasB[2] + rasB[1]*rasB[2]), data_(new DataType[size_]), rasA_(rasA), rasB_(rasB) {
      zero();
    }
    ASD_RAS_ActiveRotationMatrix(const ASD_RAS_ActiveRotationMatrix& o) : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), nactA_(o.nactA_), nactB_(o.nactB_),  size_(o.size_), data_(new DataType[o.size_]), rasA_(o.rasA_), rasB_(o.rasB_) {
      *this = o;
    }
    ASD_RAS_ActiveRotationMatrix(std::shared_ptr<const ASD_RAS_ActiveRotationMatrix> o)
      : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), nactA_(o->nactA_), nactB_(o->nactB_), size_(o->size_), data_(new DataType[o->size_]), rasA_(o->rasA_), rasB_(o->rasB_) {
      *this = *o;
    }
    //Matrix to ASDRotationMatrix conversion
    ASD_RAS_ActiveRotationMatrix(std::shared_ptr<const Matrix_base<DataType>> o, const int iclos, const int iact, const int ivirt, const int iactA, const int iactB, std::array<int,3> rasA, std::array<int,3> rasB)
      : nclosed_(iclos), nact_(iact), nvirt_(ivirt), nactA_(iactA), nactB_(iactB), size_(iactA*iactB
        + rasA[0]*rasA[1] + rasA[0]*rasA[2] + rasA[1]*rasA[2] 
        + rasB[0]*rasB[1] + rasB[0]*rasB[2] + rasB[1]*rasB[2]), data_(new DataType[size_]), rasA_(rasA), rasB_(rasB) {

      assert(iact == iactA+iactB);
      for (int i = 0; i != nactB_; ++i) {
        for (int j = 0; j != nactA_; ++j) {
          ele_aa(j, i) = o->element(j+nclosed_,i+nclosed_+nactA_);
        }
      }
      //TODO: template aaXX 12, 13, 23 with <0> for A <1> for B
      //A
      for (int i = 0; i != rasA_[1]; ++i) 
        for (int j = 0; j != rasA_[0]; ++ j) {
          ele_aa12A(j, i) = o->element(j+nclosed_, i+nclosed_+rasA_[0]);
        }
      for (int i = 0; i != rasA_[2]; ++i) 
        for (int j = 0; j != rasA_[0]; ++ j) {
          ele_aa13A(j, i) = o->element(j+nclosed_, i+nclosed_+rasA_[0]+rasA_[1]);
        }
      for (int i = 0; i != rasA_[2]; ++i) 
        for (int j = 0; j != rasA_[1]; ++ j) {
          ele_aa23A(j, i) = o->element(j+nclosed_+rasA_[0], i+nclosed_+rasA_[0]+rasA_[1]);
        }
      //B
      for (int i = 0; i != rasB_[1]; ++i) 
        for (int j = 0; j != rasB_[0]; ++ j) {
          ele_aa12A(j, i) = o->element(j+nclosed_+nactA_, i+nclosed_+nactA_+rasB_[0]);
        }
      for (int i = 0; i != rasB_[2]; ++i) 
        for (int j = 0; j != rasB_[0]; ++ j) {
          ele_aa13A(j, i) = o->element(j+nclosed_+nactA_, i+nclosed_+nactA_+rasB_[0]+rasB_[1]);
        }
      for (int i = 0; i != rasB_[2]; ++i) 
        for (int j = 0; j != rasB_[1]; ++ j) {
          ele_aa23A(j, i) = o->element(j+nclosed_+nactA_+rasB_[0], i+nclosed_+nactA_+rasB_[0]+rasB_[1]);
        }
    }

    std::shared_ptr<ASD_RAS_ActiveRotationMatrix<DataType>> clone() const {
      return std::make_shared<ASD_RAS_ActiveRotationMatrix<DataType>>(nclosed_, nact_, nvirt_, nactA_, nactB_, rasA_, rasB_);
    }
    std::shared_ptr<ASD_RAS_ActiveRotationMatrix<DataType>> copy() const {
      return std::make_shared<ASD_RAS_ActiveRotationMatrix<DataType>>(*this);
    }

    // overloaded operators
    ASD_RAS_ActiveRotationMatrix<DataType> operator+(const ASD_RAS_ActiveRotationMatrix<DataType>& o) const { ASD_RAS_ActiveRotationMatrix<DataType> out(*this); out.ax_plus_y(1.0, o); return out; }
    ASD_RAS_ActiveRotationMatrix<DataType> operator-(const ASD_RAS_ActiveRotationMatrix<DataType>& o) const { ASD_RAS_ActiveRotationMatrix<DataType> out(*this); out.ax_plus_y(-1.0, o); return out; }
    ASD_RAS_ActiveRotationMatrix<DataType>& operator+=(const ASD_RAS_ActiveRotationMatrix<DataType>& o) { ax_plus_y(1.0, o); return *this; }
    ASD_RAS_ActiveRotationMatrix<DataType>& operator-=(const ASD_RAS_ActiveRotationMatrix<DataType>& o) { ax_plus_y(-1.0, o); return *this; }
    ASD_RAS_ActiveRotationMatrix<DataType>& operator*=(const DataType a) { scale(a); return *this; }
    ASD_RAS_ActiveRotationMatrix<DataType>& operator/=(const ASD_RAS_ActiveRotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)/= o.data(i); return *this; }
    ASD_RAS_ActiveRotationMatrix<DataType>& operator*=(const ASD_RAS_ActiveRotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)*= o.data(i); return *this; }
    ASD_RAS_ActiveRotationMatrix<DataType> operator/(const ASD_RAS_ActiveRotationMatrix<DataType>& o) const { ASD_RAS_ActiveRotationMatrix<DataType> out(*this); return out /= o; }
    ASD_RAS_ActiveRotationMatrix<DataType> operator*(const ASD_RAS_ActiveRotationMatrix<DataType>& o) const { ASD_RAS_ActiveRotationMatrix<DataType> out(*this); return out *= o; }
    ASD_RAS_ActiveRotationMatrix<DataType>& operator=(const ASD_RAS_ActiveRotationMatrix<DataType>& o) { std::copy_n(o.data(), size(), data());  return *this; }

    // size of the file
    int size() const { return size_; }
    // zero out
    void zero() { fill(DataType(0.0)); }
    void fill(const DataType& a) { std::fill_n(data(), size_, a); }
    // returns dot product
    DataType dot_product(const ASD_RAS_ActiveRotationMatrix<DataType>& o) const { return blas::dot_product(data(), size_, o.data()); }
    DataType dot_product(const std::shared_ptr<const ASD_RAS_ActiveRotationMatrix<DataType>> o) const { return dot_product(*o); }
    // scale function
    void scale(const DataType& a) { std::for_each(data(), data()+size_, [&a](DataType& p) { p *= a; }); }
    // returns norm of the vector
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    // returns a root mean square
    double rms() const { return norm() / std::sqrt(static_cast<double>(size())); }
    // daxpy added to self
    void ax_plus_y(const DataType& a, const ASD_RAS_ActiveRotationMatrix& o) { std::transform(o.data(), o.data()+size_, data(), data(), [&a](DataType p, DataType q) { return p*a+q; }); }
    void ax_plus_y(const DataType& a, const std::shared_ptr<const ASD_RAS_ActiveRotationMatrix> o) { ax_plus_y(a, *o); }

    // orthogonalize to the liset of ASD_RAS_ActiveRotationMatrix's
    double orthog(std::list<std::shared_ptr<const ASD_RAS_ActiveRotationMatrix<DataType>>> c) {
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
    DataType& ele_aa(const int ia, const int ib) { return data_[ia + ib*nactA_]; }
    // RAS A
    DataType* ptr_aa12A() { return data() + nactA_*nactB_; }
    DataType& ele_aa12A(const int ia, const int ib) { return data_[nactA_*nactB_ 
                                                                  + ia + ib*rasA_[0]]; }
    DataType* ptr_aa13A() { return data() + nactA_*nactB_ + rasA_[0]*rasA_[1]; }
    DataType& ele_aa13A(const int ia, const int ib) { return data_[nactA_*nactB_ 
                                                                  + rasA_[0]*rasA_[1] + ia + ib*rasA_[0]]; }
    DataType* ptr_aa23A() { return data() + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2]; }
    DataType& ele_aa23A(const int ia, const int ib) { return data_[nactA_*nactB_ 
                                                                  + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + ia + ib*rasA_[1]]; }
    // RASB
    DataType* ptr_aa12B() { return data() + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]; }
    DataType& ele_aa12B(const int ia, const int ib) { return data_[nactA_*nactB_ 
                                                                  + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] 
                                                                  + ia + ib*rasB_[0]]; }
    DataType* ptr_aa13B() { return data() + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]
                                                                                                      + rasB_[0]*rasB_[1]; }
    DataType& ele_aa13B(const int ia, const int ib) { return data_[nactA_*nactB_ 
                                                                  + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] 
                                                                  + rasB_[0]*rasB_[1] + ia + ib*rasB_[0]]; }
    DataType* ptr_aa23B() { return data() + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]
                                                                                                      + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2]; }
    DataType& ele_aa23B(const int ia, const int ib) { return data_[nactA_*nactB_ 
                                                                  + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] 
                                                                  + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2] + ia + ib*rasB_[1]]; }


    // const references and pointers
    const DataType& ele_aa(const int ia, const int ib) const { return data_[ia + ib*nactA_]; }
    const DataType& ele_aa12A(const int ia, const int ib) const { return data_[nactA_*nactB_ + ia + ib*rasA_[0]]; }
    const DataType& ele_aa13A(const int ia, const int ib) const { return data_[nactA_*nactB_ + rasA_[0]*rasA_[1] + ia + ib*rasA_[0]]; }
    const DataType& ele_aa23A(const int ia, const int ib) const { return data_[nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + ia + ib*rasA_[1]]; }
    const DataType& ele_aa12B(const int ia, const int ib) const { return data_[nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + ia + ib*rasB_[0]]; }
    const DataType& ele_aa13B(const int ia, const int ib) const { return data_[nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1] + ia + ib*rasB_[0]]; }
    const DataType& ele_aa23B(const int ia, const int ib) const { return data_[nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2] + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2] + ia + ib*rasB_[1]]; }

    const DataType* ptr_aa() const { return data(); }
    const DataType* ptr_aa12A() const { return data()  + nactA_*nactB_; }
    const DataType* ptr_aa13A() const { return data()  + nactA_*nactB_ + rasA_[0]*rasA_[1]; }
    const DataType* ptr_aa23A() const { return data()  + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2]; }
    const DataType* ptr_aa12B() const { return data()  + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]; }
    const DataType* ptr_aa13B() const { return data()  + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]
                                                                                                                  + rasB_[0]*rasB_[1]; }
    const DataType* ptr_aa23B() const { return data()  + nactA_*nactB_ + rasA_[0]*rasA_[1] + rasA_[0]*rasA_[2] + rasA_[1]*rasA_[2]
                                                                                                                  + rasB_[0]*rasB_[1] + rasB_[0]*rasB_[2]; }

    // unpack to Matrix
    template<class MatType>
    std::shared_ptr<MatType> unpack(const DataType a = 0.0) const {
      const int nbasis = nclosed_ + nact_ + nvirt_;
      auto out = std::make_shared<MatType>(nbasis, nbasis);
      std::fill_n(out->data(), out->size(), a);

      //Active-Active 
      for (int i = 0; i != nactB_; ++i) {
        for (int j = 0; j != nactA_; ++j) {
        //out->element(j+nclosed_, i+nclosed_+nactA_) = ele_aa(j, i);
          out->element(i+nclosed_+nactA_,j+nclosed_) = -ele_aa(j,i);
        //out->element(i+nclosed_+nactA_,j+nclosed_) = ele_aa(j,i);
        }
      }
      //RAS A: follows Active-active parity
      for (int i = 0; i != rasA_[1]; ++i)
        for (int j = 0; j != rasA_[0]; ++j) {
          out->element(i+nclosed_+rasA_[0],j+nclosed_) = -ele_aa12A(j,i);
        }
      for (int i = 0; i != rasA_[2]; ++i) 
        for (int j = 0; j != rasA_[0]; ++j) {
          out->element(i+nclosed_+rasA_[0]+rasA_[1],j+nclosed_) = -ele_aa13A(j,i);
        }
      for (int i = 0; i != rasA_[2]; ++i) 
        for (int j = 0; j != rasA_[1]; ++j) {
          out->element(i+nclosed_+rasA_[0]+rasA_[1],j+nclosed_+rasA_[0]) = -ele_aa23A(j,i);
        }
      //RAS B
      for (int i = 0; i != rasB_[1]; ++i) 
        for (int j = 0; j != rasB_[0]; ++j) {
          out->element(i+nclosed_+nactA_+rasB_[0],j+nclosed_+nactA_) = -ele_aa12B(j,i);
        }
      for (int i = 0; i != rasB_[2]; ++i) 
        for (int j = 0; j != rasB_[0]; ++j) {
          out->element(i+nclosed_+nactA_+rasB_[0]+rasB_[1],j+nclosed_+nactA_) = -ele_aa13B(j,i);
        }
      for (int i = 0; i != rasB_[2]; ++i) 
        for (int j = 0; j != rasB_[1]; ++j) {
          out->element(i+nclosed_+nactA_+rasB_[0]+rasB_[1],j+nclosed_+nactA_+rasB_[0]) = -ele_aa23B(j,i);
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
        std::cout << " printing active(A)-active(B) block" << std::endl;
        for (int i = 0; i != nactB_; ++i) {
          for (int j = 0; j != nactA_; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa(j,i);
          }
          std::cout << std::endl;
        }
      }
      if (nactA_ && nactB_) {
        std::cout << " printing RAS(A) blocks" << std::endl;
        std::cout << "          RAS1-RAS2 blocks" << std::endl;
        for (int i = 0; i != rasA_[1]; ++i) 
          for (int j = 0; j != rasA_[0]; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa12A(j,i);
          }
        std::cout << std::endl;
        std::cout << "          RAS1-RAS3 blocks" << std::endl;
        for (int i = 0; i != rasA_[2]; ++i) 
          for (int j = 0; j != rasA_[0]; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa13A(j,i);
          }
        std::cout << std::endl;
        std::cout << "          RAS2-RAS3 blocks" << std::endl;
        for (int i = 0; i != rasA_[2]; ++i) 
          for (int j = 0; j != rasA_[1]; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa23A(j,i);
          }
        std::cout << std::endl;
        std::cout << " printing RAS(B) blocks" << std::endl;
        std::cout << "          RAS1-RAS2 blocks" << std::endl;
        for (int i = 0; i != rasB_[1]; ++i) 
          for (int j = 0; j != rasB_[0]; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa12B(j,i);
          }
        std::cout << std::endl;
        std::cout << "          RAS1-RAS3 blocks" << std::endl;
        for (int i = 0; i != rasB_[2]; ++i) 
          for (int j = 0; j != rasB_[0]; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa13B(j,i);
          }
        std::cout << std::endl;
        std::cout << "          RAS2-RAS3 blocks" << std::endl;
        for (int i = 0; i != rasB_[2]; ++i) 
          for (int j = 0; j != rasB_[1]; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa23B(j,i);
          }
        std::cout << std::endl;
      }
    }
};

using ASD_RAS_ActiveRotFile = ASD_RAS_ActiveRotationMatrix<double>;

}

#endif
