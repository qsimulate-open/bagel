//
// BAGEL - Parallel electron correlation program.
// Filename: rotationmatrix.h
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

//TODO

#ifndef __BAGEL_ASDSCF_ROTATIONMATRIX_H
#define __BAGEL_ASDSCF_ROTATIONMATRIX_H

#include <src/math/matrix.h>

namespace bagel {

class RotationMatrix {
  protected:
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    const int nactA_;
    const int nactB_;
    const int size_;
    std::unique_ptr<double[]> data_;

  public:
    RotationMatrix(const int iclos, const int iact, const int ivirt, const int iactA, const int iactB)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), nactA_(iactA), nactB_(iactB), size_(iclos*iact+iclos*ivirt+iact*ivirt+iactA*iactB), data_(new double[size_]) {
      zero();
    }

    RotationMatrix(const RotationMatrix& o) : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), nactA_(o.nactA_), nactB_(o.nactB_),  size_(o.size_), data_(new double[o.size_]) {
      *this = o;
    }

    RotationMatrix(std::shared_ptr<const RotationMatrix> o)
      : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), nactA_(o->nactA_), nactB_(o->nactB_), size_(o->size_), data_(new double[o->size_]) {
      *this = *o;
    }

    //Matrix to ASDRotationMatrix conversion
    RotationMatrix(std::shared_ptr<const Matrix_base<double>> o, const int iclos, const int iact, const int ivirt, const int iactA, const int iactB)
      : nclosed_(iclos), nact_(iact), nvirt_(ivirt), nactA_(iactA), nactB_(iactB), size_(iclos*iact+iclos*ivirt+iact*ivirt+iactA*iactB), data_(new double[size_]) {
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
      for (int i = 0; i != nactB_; ++i) {
        for (int j = 0; j != nactA_; ++j) {
          ele_aa(j, i) = o->element(j+nclosed_,i+nclosed_+nactA_);
        }
      }
    }

    std::shared_ptr<RotationMatrix> clone() const {
      return std::make_shared<RotationMatrix>(nclosed_, nact_, nvirt_, nactA_, nactB_);
    }

    std::shared_ptr<RotationMatrix> copy() const {
      return std::make_shared<RotationMatrix>(*this);
    }

    // overloaded operators
    RotationMatrix  operator+ (const RotationMatrix& o) const { RotationMatrix out(*this); out.ax_plus_y(1.0, o); return out; }
    RotationMatrix  operator- (const RotationMatrix& o) const { RotationMatrix out(*this); out.ax_plus_y(-1.0, o); return out; }
    RotationMatrix& operator+=(const RotationMatrix& o) { ax_plus_y(1.0, o); return *this; }
    RotationMatrix& operator-=(const RotationMatrix& o) { ax_plus_y(-1.0, o); return *this; }
    RotationMatrix& operator*=(const double a) { scale(a); return *this; }
    RotationMatrix& operator/=(const RotationMatrix& o) { for (int i = 0; i != size(); ++i) data(i)/= o.data(i); return *this; }
    RotationMatrix& operator*=(const RotationMatrix& o) { for (int i = 0; i != size(); ++i) data(i)*= o.data(i); return *this; }
    RotationMatrix  operator/ (const RotationMatrix& o) const { RotationMatrix out(*this); return out /= o; }
    RotationMatrix  operator* (const RotationMatrix& o) const { RotationMatrix out(*this); return out *= o; }
    RotationMatrix& operator= (const RotationMatrix& o) { std::copy_n(o.data(), size(), data());  return *this; }

    // size of the file
    int size() const { return size_; }
    // zero out
    void zero() { fill(double(0.0)); }
    void fill(const double& a) { std::fill_n(data(), size_, a); }
    // returns dot product
    double dot_product(const RotationMatrix& o) const { return blas::dot_product(data(), size_, o.data()); }
    double dot_product(const std::shared_ptr<const RotationMatrix> o) const { return dot_product(*o); }
    // scale function
    void scale(const double& a) { std::for_each(data(), data()+size_, [&a](double& p) { p *= a; }); }
    // returns norm of the vector
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    // returns a root mean square
    double rms() const { return norm() / std::sqrt(static_cast<double>(size())); }
    // daxpy added to self
    void ax_plus_y(const double& a, const RotationMatrix& o) { std::transform(o.data(), o.data()+size_, data(), data(), [&a](double p, double q) { return p*a+q; }); }
    void ax_plus_y(const double& a, const std::shared_ptr<const RotationMatrix> o) { ax_plus_y(a, *o); }

    // orthogonalize to the liset of RotationMatrix's
    double orthog(std::list<std::shared_ptr<const RotationMatrix>> c) {
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
    double* data() { return data_.get(); }
    const double* data() const { return data_.get(); }
    double& data(const size_t i) { return data_[i]; }
    const double& data(const size_t i) const { return data_[i]; }
    // return data_
    double* begin() { return data(); }
    // return data_
    double* end() { return data()+size_; }

    // ROTATION MATRIX
    //   C   A   V             A(A) A(B)   Stored in order of:
    // C .  ca  cv  where  A(A)  .  [aa]   ca -> va -> vc -> aa
    // A . [aa]  .         A(B)  .    .
    // V .  va   .

    // closed-active block. closed runs first
    double* ptr_ca() { return data(); }                     //row  col
    double& ele_ca(const int ic, const int ia) { return data_[ic + ia*nclosed_]; }
    // active-virtual block. virtual runs first
    double* ptr_va() { return data() + nclosed_*nact_; }                     //row  col
    double& ele_va(const int iv, const int ia) { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    // closed-virtual block. virtual runs first
    double* ptr_vc() { return data() + (nclosed_+nvirt_)*nact_; }                     //row  col
    double& ele_vc(const int iv, const int ic) { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    // active-active block. A index runs first
    double* ptr_aa() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_; }                     //row  col
    double& ele_aa(const int ia, const int ib) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ia + ib*nactA_]; }

    // const references and pointers
    const double& ele_ca(const int ic, const int ia) const { return data_[ic + ia*nclosed_]; }
    const double& ele_va(const int iv, const int ia) const { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    const double& ele_vc(const int iv, const int ic) const { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    const double& ele_aa(const int ia, const int ib) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ia + ib*nactA_]; }

    // pointers to the beginning of the given block
    const double* ptr_ca() const { return data(); }
    const double* ptr_va() const { return data() + nclosed_*nact_; }
    const double* ptr_vc() const { return data() + (nclosed_+nvirt_)*nact_; }
    const double* ptr_aa() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_; }

    // unpack to Matrix:
    // fill the upper-triangular part then anti-symmetrize
    std::shared_ptr<Matrix> unpack(const double a = 0.0) const {
      const int nocc = nclosed_ + nact_;
      const int nbasis = nclosed_ + nact_ + nvirt_;
      auto out = std::make_shared<Matrix>(nbasis, nbasis);
      std::fill_n(out->data(), out->size(), a);
      //Virtual-Active
      for (int j = 0; j != nact_; ++j) {
        for (int i = 0; i != nvirt_; ++i) {
          out->element(i+nocc,j+nclosed_) = ele_va(i,j);
        } //i
      } //j

      //Active-Closed
      for (int i = 0; i != nact_; ++i) {
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
      for (int i = 0; i != nactB_; ++i) {
        for (int j = 0; j != nactA_; ++j) {
        //out->element(j+nclosed_, i+nclosed_+nactA_) = ele_aa(j, i);
          out->element(i+nclosed_+nactA_,j+nclosed_) = ele_aa(j, i);
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
        std::cout << " printing active(A)-active(B) block" << std::endl;
        for (int i = 0; i != nactB_; ++i) {
          for (int j = 0; j != nactA_; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa(j,i);
          }
          std::cout << std::endl;
        }
      }

    }
};


}

#endif
