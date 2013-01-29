//
// BAGEL - Parallel electron correlation program.
// Filename: rdm.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __NEWINT_WFN_RDM_H
#define __NEWINT_WFN_RDM_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <type_traits>
#include <src/wfn/geometry.h>
#include <src/util/f77.h>
#include <src/util/matrix.h>

namespace bagel {

class RDM_base {
  protected:
    std::unique_ptr<double[]> data_;
    const int norb_;
    size_t dim_;
    int rank_;

  public:
    RDM_base(const int n, const int rank);
    RDM_base(const RDM_base& o);

    double* first() { return data(); }
    double* data() { return data_.get(); }
    const double* data() const { return data_.get(); }

    void zero() { std::fill(data(), data()+dim_*dim_, 0.0); }
    void daxpy(const double a, const RDM_base& o) { daxpy_(dim_*dim_, a, o.data(), 1, data(), 1); }
    void daxpy(const double a, const std::shared_ptr<RDM_base>& o) { this->daxpy(a, *o); }
    void scale(const double a) { dscal_(dim_*dim_, a, data(), 1); }
    size_t size() const { return dim_*dim_; }

    int norb() const { return norb_; }

};


template <int rank>
class RDM : public RDM_base {
  protected:
    // T should be able to be multiplied by norb_
    template<int i, typename T, typename ...args>
    size_t address_(const T& head, const args&... tail) const {
      static_assert(i >= 0 && std::is_integral<T>::value, "address_ called with a wrong template variable");
      T out = head;
      for (int j = 0; j != i; ++j) out *= norb_; 
      return out + address_<i+1>(tail...);
    }
    template<int i, typename T>
    size_t address_(const T& head) const {
      static_assert(i+1 == rank*2 && std::is_integral<T>::value, "address_(const T&) called with a wrong template variable");
      T out = head;
      for (int j = 0; j != i; ++j) out *= norb_; 
      return out;
    }

  public:
    RDM(const int n) : RDM_base(n, rank) { };
    RDM(const RDM& o) : RDM_base(o) {};
    ~RDM() {  };

    std::shared_ptr<RDM<rank>> clone() const { return std::shared_ptr<RDM<rank>>(new RDM<rank>(norb_)); };

    template<typename ...args>
    double& element(const args&... index) { return data_[address_<0>(index...)]; }

    template<typename ...args>
    const double& element(const args&... index) const { return data_[address_<0>(index...)]; }


    // returns if this is natural orbitals - only for rank 1
    bool natural_orbitals() const {
      if (rank != 1) throw std::logic_error("RDM::natural_orbitals() is implemented only for rank 1");
      const double a = ddot_(norb_*norb_, data_, 1, data_, 1);
      const double b = ddot_(norb_, data_, norb_+1, data_, norb_+1);
      return ::fabs(a-b) < 1.0e-12;
    }


    std::shared_ptr<Matrix> rdm1_mat(std::shared_ptr<const Geometry> g, const int nclosed, const bool all = true) const {
      static_assert(rank == 1, "RDM::rdm1_mat is only implemented for rank == 1");
      std::shared_ptr<Matrix> out(new Matrix(nclosed+norb_, nclosed+norb_));
      if (all)
        for (int i = 0; i != nclosed; ++i) out->element(i,i) = 2.0;
      for (int i = 0; i != norb_; ++i)
        for (int j = 0; j != norb_; ++j)
          out->element(j+nclosed, i+nclosed) = element(j,i);
      return out;
    }

    std::pair<std::shared_ptr<Matrix>, std::vector<double>> generate_natural_orbitals() const {
      static_assert(rank == 1, "RDM::generate_natural_orbitals is only implemented for rank == 1");
      std::shared_ptr<Matrix> buf(new Matrix(dim_,dim_,true));
      buf->add_diag(2.0);
      daxpy_(dim_*dim_, -1.0, data(), 1, buf->data(), 1);

      std::vector<double> vec(dim_);
      buf->diagonalize(&vec[0]);

      for (auto& i : vec) i = 2.0-i;
      return std::make_pair(buf, vec);
    }

    void transform(const std::shared_ptr<Matrix>& coeff) {
      const double* start = coeff->data();
      std::unique_ptr<double[]> buf(new double[dim_*dim_]);
      if (rank == 1) {
        dgemm_("N", "N", dim_, dim_, dim_, 1.0, data(), dim_, start, dim_, 0.0, buf.get(), dim_);
        dgemm_("T", "N", dim_, dim_, dim_, 1.0, start, dim_, buf.get(), dim_, 0.0, data(), dim_);
      } else if (rank == 2) {
        // first half transformation
        dgemm_("N", "N", dim_*norb_, norb_, norb_, 1.0, data(), dim_*norb_, start, norb_, 0.0, buf.get(), dim_*norb_);
        for (int i = 0; i != norb_; ++i)
          dgemm_("N", "N", dim_, norb_, norb_, 1.0, buf.get()+i*dim_*norb_, dim_, start, norb_, 0.0, data()+i*dim_*norb_, dim_);
        // then tranpose
        mytranspose_(data(), dim_, dim_, buf.get());
        // and do it again
        dgemm_("N", "N", dim_*norb_, norb_, norb_, 1.0, buf.get(), dim_*norb_, start, norb_, 0.0, data(), dim_*norb_);
        for (int i = 0; i != norb_; ++i)
          dgemm_("N", "N", dim_, norb_, norb_, 1.0, data()+i*dim_*norb_, dim_, start, norb_, 0.0, buf.get()+i*dim_*norb_, dim_);
        // to make sure for non-symmetric density matrices (and anyway this should be cheap).
        mytranspose_(buf.get(), dim_, dim_, data());
      } else {
        assert(false);
      }
    }

    std::vector<double> diag() const {
      std::vector<double> out(dim_);
      for (int i = 0; i != dim_; ++i) out[i] = element(i,i);
      return out;
    }


    // What is the best way to implement this for general rank??
    void print(const double thresh = 1.0e-3) const {
      static_assert(rank <= 3, "RDM::print is so far only implemented for RDM1 and 2");
      const double* ptr = data_.get();
      if (rank == 1) {
        for (int i = 0; i != norb_; ++i) {
          for (int j = 0; j != norb_; ++j, ++ptr)
            std::cout << std::setw(12) << std::setprecision(7) << *ptr;
          std::cout << std::endl;
        }
      } else if (rank == 2) {
        for (int i = 0; i != norb_; ++i) {
          for (int j = 0; j != norb_; ++j) {
            for (int k = 0; k != norb_; ++k) {
              for (int l = 0; l != norb_; ++l, ++ptr) {
                if (std::abs(*ptr) > thresh) std::cout << std::setw(3) << l << std::setw(3)
                      << k << std::setw(3) << j << std::setw(3) << i
                      << std::setw(12) << std::setprecision(7) << *ptr << std::endl;
        } } } }
      } else if (rank == 3) {
        for (int i = 0; i != norb_; ++i) {
          for (int j = 0; j != norb_; ++j) {
            for (int k = 0; k != norb_; ++k) {
              for (int l = 0; l != norb_; ++l) {
              for (int m = 0; m != norb_; ++m) {
              for (int n = 0; n != norb_; ++n, ++ptr) {
                if (std::abs(*ptr) > thresh) std::cout << std::setw(3) << n << std::setw(3) << m << std::setw(3) << l << std::setw(3)
                      << k << std::setw(3) << j << std::setw(3) << i
                      << std::setw(12) << std::setprecision(7) << *ptr << std::endl;
        } } } } } }
      }
    }
};

}

#endif
