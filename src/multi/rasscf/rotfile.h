//
// BAGEL - Parallel electron correlation program.
// Filename: multi/rasscf/rotfile.h
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


#ifndef __BAGEL_RASSCF_ROTFILE_H
#define __BAGEL_RASSCF_ROTFILE_H

#include <src/util/math/matrix.h>

namespace bagel {

template<typename DataType>
class RASRotationMatrix {
  public:
    using data_type = DataType;

  protected:
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    std::array<int,3> ras_;
    const int size_;
    std::unique_ptr<DataType[]> data_;

  public:
    RASRotationMatrix(const int iclos, const int iact, const int ivirt, std::array<int,3> rsp)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt), ras_(rsp), 
       size_(iclos*iact+iclos*ivirt+iact*ivirt + rsp[0]*rsp[1] + rsp[0]*rsp[2] + rsp[1]*rsp[2] ), 
       data_(new DataType[size_]) {
      zero();
    }
    RASRotationMatrix(const RASRotationMatrix& o) : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), ras_(o.ras_), size_(o.size_), data_(new DataType[o.size_]) {
      *this = o;
    }
    RASRotationMatrix(std::shared_ptr<const RASRotationMatrix> o)
      : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), ras_(o->ras_), size_(o->size_), data_(new DataType[o->size_]) {
      *this = *o;
    }
    RASRotationMatrix(std::shared_ptr<const Matrix_base<DataType>> o, const int iclos, const int iact, const int ivirt, std::array<int,3> rsp)
      : nclosed_(iclos), nact_(iact), nvirt_(ivirt), ras_(rsp), 
        size_(iclos*iact+iclos*ivirt+iact*ivirt + rsp[0]*rsp[1] + rsp[0]*rsp[2] + rsp[1]*rsp[2] ), 
        data_(new DataType[size_]) {
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
      //RAS
      for (int i = 0; i != ras_[1]; ++i)
        for (int j = 0; j != ras_[0]; ++ j) {
          ele_aa12(j, i) = o->element(j+nclosed_, i+nclosed_+ras_[0]);
        }
      for (int i = 0; i != ras_[2]; ++i)
        for (int j = 0; j != ras_[0]; ++ j) {
          ele_aa13(j, i) = o->element(j+nclosed_, i+nclosed_+ras_[0]+ras_[1]);
        }
      for (int i = 0; i != ras_[2]; ++i)
        for (int j = 0; j != ras_[1]; ++ j) {
          ele_aa23(j, i) = o->element(j+nclosed_+ras_[0], i+nclosed_+ras_[0]+ras_[1]);
        }
    }

    std::shared_ptr<RASRotationMatrix<DataType>> clone() const {
      return std::make_shared<RASRotationMatrix<DataType>>(nclosed_, nact_, nvirt_, ras_);
    }
    std::shared_ptr<RASRotationMatrix<DataType>> copy() const {
      return std::make_shared<RASRotationMatrix<DataType>>(*this);
    }

    // overloaded operators
    RASRotationMatrix<DataType> operator+(const RASRotationMatrix<DataType>& o) const { RASRotationMatrix<DataType> out(*this); out.ax_plus_y(1.0, o); return out; }
    RASRotationMatrix<DataType> operator-(const RASRotationMatrix<DataType>& o) const { RASRotationMatrix<DataType> out(*this); out.ax_plus_y(-1.0, o); return out; }
    RASRotationMatrix<DataType>& operator+=(const RASRotationMatrix<DataType>& o) { ax_plus_y(1.0, o); return *this; }
    RASRotationMatrix<DataType>& operator-=(const RASRotationMatrix<DataType>& o) { ax_plus_y(-1.0, o); return *this; }
    RASRotationMatrix<DataType>& operator*=(const DataType a) { scale(a); return *this; }
    RASRotationMatrix<DataType>& operator/=(const RASRotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)/= o.data(i); return *this; }
    RASRotationMatrix<DataType>& operator*=(const RASRotationMatrix<DataType>& o) { for (int i = 0; i != size(); ++i) data(i)*= o.data(i); return *this; }
    RASRotationMatrix<DataType> operator/(const RASRotationMatrix<DataType>& o) const { RASRotationMatrix<DataType> out(*this); return out /= o; }
    RASRotationMatrix<DataType> operator*(const RASRotationMatrix<DataType>& o) const { RASRotationMatrix<DataType> out(*this); return out *= o; }
    RASRotationMatrix<DataType>& operator=(const RASRotationMatrix<DataType>& o) { std::copy_n(o.data(), size(), data());  return *this; }

    // size of the file
    int size() const { return size_; }
    // zero out
    void zero() { fill(DataType(0.0)); }
    void fill(const DataType& a) { std::fill_n(data(), size_, a); }
    // returns dot product
    DataType dot_product(const RASRotationMatrix<DataType>& o) const { return blas::dot_product(data(), size_, o.data()); }
    DataType dot_product(const std::shared_ptr<const RASRotationMatrix<DataType>> o) const { return dot_product(*o); }
    // scale function
    void scale(const DataType& a) { std::for_each(data(), data()+size_, [&a](DataType& p) { p *= a; }); }
    // returns norm of the vector
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    // returns a root mean square
    double rms() const { return norm() / std::sqrt(static_cast<double>(size())); }
    // daxpy added to self
    void ax_plus_y(const DataType& a, const RASRotationMatrix& o) { std::transform(o.data(), o.data()+size_, data(), data(), [&a](DataType p, DataType q) { return p*a+q; }); }
    void ax_plus_y(const DataType& a, const std::shared_ptr<const RASRotationMatrix> o) { ax_plus_y(a, *o); }

    // orthogonalize to the liset of RASRotationMatrix's
    double orthog(std::list<std::shared_ptr<const RASRotationMatrix<DataType>>> c) {
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

// RAS 
    DataType& ele_aa12(const int ia, const int ib) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ia + ib*ras_[0]]; }
    DataType& ele_aa13(const int ia, const int ib) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ras_[0]*ras_[1] + ia + ib*ras_[0]]; }
    DataType& ele_aa23(const int ia, const int ib) { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ras_[0]*ras_[1] + ras_[0]*ras_[2] + ia + ib*ras_[1]]; }
    DataType* ptr_aa12() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_; }
    DataType* ptr_aa13() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ras_[0]*ras_[1]; }
    DataType* ptr_aa23() { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ras_[0]*ras_[1] + ras_[0]*ras_[2]; }

    const DataType& ele_aa12(const int ia, const int ib) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ia + ib*ras_[0]]; }
    const DataType& ele_aa13(const int ia, const int ib) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ras_[0]*ras_[1] + ia + ib*ras_[0]]; }
    const DataType& ele_aa23(const int ia, const int ib) const { return data_[(nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ras_[0]*ras_[1] + ras_[0]*ras_[2] + ia + ib*ras_[1]]; }
    const DataType* ptr_aa12() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_; }
    const DataType* ptr_aa13() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ras_[0]*ras_[1]; }
    const DataType* ptr_aa23() const { return data() + (nclosed_+nvirt_)*nact_ + nclosed_*nvirt_ + ras_[0]*ras_[1] + ras_[0]*ras_[2]; }

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
      //RAS 
      for (int i = 0; i != ras_[1]; ++i)
        for (int j = 0; j != ras_[0]; ++j) {
          out->element(i+nclosed_+ras_[0],j+nclosed_) = -ele_aa12(j,i);
        }
      for (int i = 0; i != ras_[2]; ++i)
        for (int j = 0; j != ras_[0]; ++j) {
          out->element(i+nclosed_+ras_[0]+ras_[1],j+nclosed_) = -ele_aa13(j,i);
        }
      for (int i = 0; i != ras_[2]; ++i)
        for (int j = 0; j != ras_[1]; ++j) {
          out->element(i+nclosed_+ras_[0]+ras_[1],j+nclosed_+ras_[0]) = -ele_aa23(j,i);
        }
      out->print("Rotation Matrix before anti-sym",nbasis);
      for (int i = 0; i != nbasis; ++i) {
        for (int j = 0; j <= i; ++j) {
          out->element(j, i) = - detail::conj(out->element(i, j));
        }
      }
      out->print("Rotation Matrix after anti-sym",nbasis);
      return out;
    }

    template<class MatType>
    std::shared_ptr<MatType> unpack_sym(const DataType a = 0.0) const {
      assert(false);
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

      if (nact_) {
        std::cout << " printing RAS blocks" << std::endl;
        std::cout << "          RAS1-RAS2 blocks" << std::endl;
        for (int i = 0; i != ras_[1]; ++i)
          for (int j = 0; j != ras_[0]; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa12(j,i);
          }
        std::cout << std::endl;
        std::cout << "          RAS1-RAS3 blocks" << std::endl;
        for (int i = 0; i != ras_[2]; ++i)
          for (int j = 0; j != ras_[0]; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa13(j,i);
          }
        std::cout << std::endl;
        std::cout << "          RAS2-RAS3 blocks" << std::endl;
        for (int i = 0; i != ras_[2]; ++i)
          for (int j = 0; j != ras_[1]; ++j) {
            std::cout << std::setw(10) << std::setprecision(6) << ele_aa23(j,i);
          }
        std::cout << std::endl;
      }

    }
};

using RASRotFile = RASRotationMatrix<double>;
//using ZRASRotFile = RASRotationMatrix<std::complex<double>>;

}

#endif
