//
// BAGEL - Parallel electron correlation program.
// Filename: rotfile.h
// Copyright (C) 2011 Toru Shiozaki
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


#ifndef __BAGEL_CASSCF_ROTFILE_H
#define __BAGEL_CASSCF_ROTFILE_H

#include <src/math/matrix.h>

namespace bagel {

template<typename DataType>
class RotationMatrix {
  protected:
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    const bool superci_;
    const int size_;
    std::unique_ptr<DataType[]> data_;

  public:
    RotationMatrix(const int iclos, const int iact, const int ivirt, const bool superci = true)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), superci_(superci), size_(iclos*iact+iclos*ivirt+iact*ivirt+(superci ? 1 : 0)), data_(new DataType[size_]) {
      zero();
    }
    RotationMatrix(const RotationMatrix& o) : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), superci_(o.superci_), size_(o.size_), data_(new DataType[o.size_]) {
      *this = o;
    }
    RotationMatrix(std::shared_ptr<const RotationMatrix> o)
      : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), superci_(o->superci_), size_(o->size_), data_(new DataType[o->size_]) {
      *this = *o;
    }
    RotationMatrix(std::shared_ptr<const Matrix_base<DataType>> o, const int iclos, const int iact, const int ivirt, const bool superci = true)
      : nclosed_(iclos), nact_(iact), nvirt_(ivirt), superci_(superci), size_(iclos*iact+iclos*ivirt+iact*ivirt+(superci ? 1 : 0)), data_(new DataType[size_]) {
      const int nocc = nclosed_ + nact_;
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
    }

    std::shared_ptr<RotationMatrix<DataType>> clone() const {
      return std::make_shared<RotationMatrix<DataType>>(nclosed_, nact_, nvirt_, superci_);
    }
    std::shared_ptr<RotationMatrix<DataType>> copy() const {
      return std::make_shared<RotationMatrix<DataType>>(*this);
    }

    // overloaded operators
    RotationMatrix<DataType> operator+(const RotationMatrix<DataType>& o) const { RotationMatrix<DataType> out(*this); out.ax_plus_y(1.0, o); return out; }
    RotationMatrix<DataType> operator-(const RotationMatrix<DataType>& o) const { RotationMatrix<DataType> out(*this); out.ax_plus_y(-1.0, o); return out; }
    RotationMatrix<DataType>& operator+=(const RotationMatrix<DataType>& o) { ax_plus_y(1.0, o); return *this; }
    RotationMatrix<DataType>& operator-=(const RotationMatrix<DataType>& o) { ax_plus_y(-1.0, o); return *this; }
    RotationMatrix<DataType>& operator*=(const DataType a) { scale(a); return *this; }
    RotationMatrix<DataType>& operator/=(const RotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)/= o.data(i); return *this; }
    RotationMatrix<DataType> operator/(const RotationMatrix<DataType>& o) const { RotationMatrix<DataType> out(*this); return out /= o; }
    RotationMatrix<DataType>& operator=(const RotationMatrix<DataType>& o) { std::copy_n(o.data(), size(), data());  return *this; }

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
    // daxpy added to self
    void ax_plus_y(const DataType& a, const RotationMatrix& o) { std::transform(o.data(), o.data()+size_, data(), data(), [&a](DataType p, DataType q) { return p*a+q; }); }
    void ax_plus_y(const DataType& a, const std::shared_ptr<const RotationMatrix> o) { ax_plus_y(a, *o); }

    // orthogonalize to the liset of RotationMatrix's
    double orthog(std::list<std::shared_ptr<const RotationMatrix<DataType>>> c) {
      for (auto iter = c.begin(); iter != c.end(); ++iter)
        this->ax_plus_y(- this->dot_product(**iter), **iter);
      const double scal = 1.0/this->norm();
      scale(scal);
      return 1.0/scal;
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

    // closed-active block. closed runs first
    DataType* ptr_ca() { return data(); }
    DataType& ele_ca(const int ic, const int ia) { return data_[ic + ia*nclosed_]; }
    // active-virtual block. virtual runs first
    DataType* ptr_va() { return data() + nclosed_*nact_; }
    DataType& ele_va(const int iv, const int ia) { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    // closed-virtual block. virtual runs first
    DataType* ptr_vc() { return data() + (nclosed_+nvirt_)*nact_; }
    DataType& ele_vc(const int iv, const int ic) { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    // reference config.
    DataType& ele_ref() { assert(superci_); return data_[size_-1]; }
    // const references
    const DataType& ele_ca(const int ic, const int ia) const { return data_[ic + ia*nclosed_]; }
    const DataType& ele_va(const int iv, const int ia) const { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    const DataType& ele_vc(const int iv, const int ic) const { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    const DataType& ele_ref() const { assert(superci_); return data_[size_-1]; }

    // unpack to Matrix
    template<class MatType>
    std::shared_ptr<MatType> unpack(const DataType a = 0.0) const {
      const int nocc = nclosed_ + nact_;
      const int nbasis = nclosed_ + nact_ + nvirt_;
      auto out = std::make_shared<MatType>(nbasis, nbasis);
      std::fill_n(out->data(), out->size(), a);

      for (int i = 0; i != nact_; ++i) {
        for (int j = 0; j != nvirt_;   ++j) {
          out->element(j+nocc, i+nclosed_) = ele_va(j, i);
        }
        for (int j = 0; j != nclosed_; ++j) {
          out->element(i+nclosed_, j) = ele_ca(j, i);
        }
      }
      for (int i = 0; i != nclosed_; ++i) {
        for (int j = 0; j != nvirt_;   ++j) {
          out->element(j+nocc, i) = ele_vc(j, i);
        }
      }
      for (int i = 0; i != nbasis; ++i) {
        for (int j = 0; j <= i; ++j) {
          out->element(j, i) = - detail::conj(out->element(i, j));
        }
      }
      return out;
    }

    template<class MatType>
    std::shared_ptr<MatType> unpack_sym(const DataType a = 0.0) const {
      const int nocc = nclosed_ + nact_;
      const int nbasis = nclosed_ + nact_ + nvirt_;
      auto out = std::make_shared<MatType>(nbasis, nbasis);
      std::fill_n(out->data(), out->size(), a);
      for (int i = 0; i != nact_; ++i) {
        for (int j = 0; j != nvirt_;   ++j) {
          out->element(j+nocc, i+nclosed_) = ele_va(j, i);
        }
        for (int j = 0; j != nclosed_; ++j) {
          out->element(i+nclosed_, j) = ele_ca(j, i);
        }
      }
      for (int i = 0; i != nclosed_; ++i) {
        for (int j = 0; j != nvirt_;   ++j) {
          out->element(j+nocc, i) = ele_vc(j, i);
        }
      }
      for (int i = 0; i != nbasis; ++i) {
        for (int j = 0; j <= i; ++j) {
          out->element(j, i) = detail::conj(out->element(i, j));
        }
      }
      return out;
    }

    // print matrix
    void print() const {
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
      if (superci_) {
        std::cout << "reference weight " << ele_ref() << std::endl;
      }
    }
};

using RotFile = RotationMatrix<double>;
using ZRotFile = RotationMatrix<std::complex<double>>;

}

#endif
