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
    RDM(const int n) : RDM_base(n, rank) { }
    RDM(const RDM& o) : RDM_base(o) { }

    std::shared_ptr<RDM<rank>> clone() const { return std::shared_ptr<RDM<rank>>(new RDM<rank>(norb_)); }

    template<typename ...args>
    double& element(const args&... index) { return data_[address_<0>(index...)]; }

    template<typename ...args>
    const double& element(const args&... index) const { return data_[address_<0>(index...)]; }


    // returns if this is natural orbitals - only for rank 1
    bool natural_orbitals() const {
      throw std::logic_error("RDM<N>::natural_orbitals() should not be called with N>1");
      return true;
    }

    std::shared_ptr<Matrix> rdm1_mat(std::shared_ptr<const Geometry> g, const int nclosed, const bool all = true) const {
      throw std::logic_error("RDM<N>::rdm1_mat() should not be called with N>1");
      return std::shared_ptr<Matrix>();
    }

    std::pair<std::shared_ptr<Matrix>, std::vector<double>> generate_natural_orbitals() const {
      throw std::logic_error("RDM<N>::generate_natural_orbitals() should not be called with N>1");
      return std::pair<std::shared_ptr<Matrix>, std::vector<double>>();
    }

    void transform(const std::shared_ptr<Matrix>& coeff) { throw std::logic_error("RDM<N>::transform() (N>3) not implemented yet"); }

    std::vector<double> diag() const {
      std::vector<double> out(dim_);
      for (int i = 0; i != dim_; ++i) out[i] = element(i,i);
      return out;
    }

    void print(const double thresh = 1.0e-3) const { throw std::logic_error("RDM<N>::print() (N>3) not implemented yet"); }
};

template<> bool RDM<1>::natural_orbitals() const;

template<> std::pair<std::shared_ptr<Matrix>, std::vector<double>> RDM<1>::generate_natural_orbitals() const;

template<> void RDM<1>::transform(const std::shared_ptr<Matrix>& coeff);
template<> void RDM<2>::transform(const std::shared_ptr<Matrix>& coeff);

template<> std::shared_ptr<Matrix> RDM<1>::rdm1_mat(std::shared_ptr<const Geometry> g, const int nclosed, const bool all) const;

template<> void RDM<1>::print(const double thresh) const;
template<> void RDM<2>::print(const double thresh) const;
template<> void RDM<3>::print(const double thresh) const;

}

#endif
