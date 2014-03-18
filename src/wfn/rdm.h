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


#ifndef __BAGEL_WFN_RDM_H
#define __BAGEL_WFN_RDM_H

#include <type_traits>
#include <src/wfn/geometry.h>
#include <src/math/matrix.h>

namespace bagel {

template<typename DataType>
class RDM_base {
  protected:
    int norb_;
    size_t dim_;
    int rank_;

    std::unique_ptr<DataType[]> data_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << norb_ << dim_ << rank_ << make_array(data(), size());
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> norb_ >> dim_ >> rank_;
      data_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      ar >> make_array(data(), size());
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }

  public:
    RDM_base() : norb_(0), dim_(0) { }

    RDM_base(const int n, const int rank) : norb_(n), rank_(rank) {
      assert(rank > 0);
      dim_ = 1lu;
      for (int i = 0; i != rank; ++i) dim_ *= n;
      data_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      std::fill_n(data(), size(), static_cast<DataType>(0.0));
    }

    RDM_base(const RDM_base<DataType>& o) : norb_(o.norb_), dim_(o.dim_), rank_(o.rank_) {
      data_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      std::copy_n(o.data(), size(), data());
    }

    RDM_base(RDM_base<DataType>&& o) : norb_(o.norb_), dim_(o.dim_), rank_(o.rank_), data_(std::move(o.data_)) {
    }

    DataType* first() { return data(); }
    DataType* data() { return data_.get(); }
    const DataType* data() const { return data_.get(); }

    void zero() { std::fill(data(), data()+dim_*dim_, static_cast<DataType>(0.0)); }
    void ax_plus_y(const DataType a, const RDM_base<DataType>& o) {
      assert(size() == o.size());
      std::transform(o.data(), o.data()+size(), data(), data(), [&a](DataType p, DataType q){ return q+a*p;});
    }
    void ax_plus_y(const DataType& a, const std::shared_ptr<RDM_base<DataType>>& o) { this->ax_plus_y(a, *o); }
    void scale(const DataType& a) { std::for_each(data(), data()+size(), [&a](DataType& p){ p *= a; }); }
    size_t size() const { return dim_*dim_; }

    int norb() const { return norb_; }
    size_t dim() const { return dim_; }

};


template <int rank, typename DataType = double>
class RDM : public RDM_base<DataType> {
  protected:
    // T should be able to be multiplied by norb_
    template<int i, typename T, typename ...args>
    size_t address_(const T& head, const args&... tail) const {
      static_assert(i >= 0 && std::is_integral<T>::value, "address_ called with a wrong template variable");
      T out = head;
      for (int j = 0; j != i; ++j) out *= this->norb_;
      return out + address_<i+1>(tail...);
    }
    template<int i, typename T>
    size_t address_(const T& head) const {
      static_assert(i+1 == rank*2 && std::is_integral<T>::value, "address_(const T&) called with a wrong template variable");
      T out = head;
      for (int j = 0; j != i; ++j) out *= this->norb_;
      return out;
    }

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<RDM_base<DataType>>(*this);
    }

  public:
    RDM() { }
    RDM(const int n) : RDM_base<DataType>(n, rank) { }
    RDM(const RDM<rank,DataType>& o) : RDM_base<DataType>(o) { }
    RDM(RDM<rank,DataType>&& o) : RDM_base<DataType>(std::move(o)) { }

    std::shared_ptr<RDM<rank,DataType>> clone() const { return std::make_shared<RDM<rank,DataType>>(this->norb_); }
    std::shared_ptr<RDM<rank,DataType>> copy() const { return std::make_shared<RDM<rank,DataType>>(*this); }

    template<typename ...args>
    DataType& element(const args&... index) { return this->data_[address_<0>(index...)]; }

    template<typename ...args>
    const DataType& element(const args&... index) const { return this->data_[address_<0>(index...)]; }

    RDM<rank,DataType>& operator+=(const RDM<rank,DataType>& o) { this->ax_plus_y(1.0, o); return *this; }
    RDM<rank,DataType>& operator-=(const RDM<rank,DataType>& o) { this->ax_plus_y(-1.0, o); return *this; }
    RDM<rank,DataType> operator+(const RDM<rank,DataType>& o) const { RDM<rank,DataType> out(*this); out.ax_plus_y(1.0, o); return out; }
    RDM<rank,DataType> operator-(const RDM<rank,DataType>& o) const { RDM<rank,DataType> out(*this); out.ax_plus_y(-1.0, o); return out; }

    // returns if this is natural orbitals - only for rank 1
    bool natural_orbitals() const {
      throw std::logic_error("RDM<N>::natural_orbitals() should not be called with N>1");
      return true;
    }

    std::shared_ptr<Matrix> rdm1_mat(const int nclosed, const bool all = true) const {
      throw std::logic_error("RDM<N>::rdm1_mat() should not be called with N>1");
      return nullptr;
    }

    std::pair<std::shared_ptr<Matrix>, std::vector<double>> generate_natural_orbitals() const {
      throw std::logic_error("RDM<N>::generate_natural_orbitals() should not be called with N>1");
      return std::pair<std::shared_ptr<Matrix>, std::vector<double>>();
    }

    void transform(const std::shared_ptr<Matrix>& coeff) { throw std::logic_error("RDM<N>::transform() (N>3) not implemented yet"); }

    std::vector<DataType> diag() const {
      std::vector<DataType> out(this->dim_);
      for (int i = 0; i != this->dim_; ++i) out[i] = element(i,i);
      return out;
    }

    void print(const double thresh = 1.0e-3) const { throw std::logic_error("RDM<N>::print() (N>3) not implemented yet"); }
};

template <int rank>
using ZRDM = RDM<rank, std::complex<double>>;

template<> bool RDM<1,double>::natural_orbitals() const;

template<> std::pair<std::shared_ptr<Matrix>, std::vector<double>> RDM<1,double>::generate_natural_orbitals() const;

template<> void RDM<1,double>::transform(const std::shared_ptr<Matrix>& coeff);
template<> void RDM<2,double>::transform(const std::shared_ptr<Matrix>& coeff);

template<> std::shared_ptr<Matrix> RDM<1,double>::rdm1_mat(const int nclosed, const bool all) const;

template<> void RDM<1,double>::print(const double thresh) const;
template<> void RDM<2,double>::print(const double thresh) const;
template<> void RDM<3,double>::print(const double thresh) const;

}

#endif
