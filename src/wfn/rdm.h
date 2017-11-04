//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rdm.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __BAGEL_WFN_RDM_H
#define __BAGEL_WFN_RDM_H

#include <type_traits>
#include <src/util/vec.h>
#include <src/util/kramers.h>
#include <src/wfn/geometry.h>

namespace bagel {


template <int rank, typename DataType = double>
class RDM : public btas::TensorN<DataType, rank*2> {
  private:
    constexpr const static int N = rank*2;
  public:
    using btas::TensorN<DataType, N>::data;
    using btas::TensorN<DataType, N>::begin;
    using btas::TensorN<DataType, N>::end;
    using btas::TensorN<DataType, N>::cbegin;
    using btas::TensorN<DataType, N>::cend;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<btas::TensorN<DataType,N>>(*this);
    }

  public:
    RDM() { }

    RDM(const int n) : btas::TensorN<DataType, N>(btas::CRange<N>(btas::Range1(n),N)) {
      zero();
    }

    RDM(const RDM<rank,DataType>& o) : btas::TensorN<DataType, N>(o) {
    }

    RDM(RDM<rank,DataType>&& o) : btas::TensorN<DataType, N>(std::move(o)) {
    }

    std::shared_ptr<RDM<rank,DataType>> clone() const { return std::make_shared<RDM<rank,DataType>>(this->norb()); }
    std::shared_ptr<RDM<rank,DataType>> copy() const { return std::make_shared<RDM<rank,DataType>>(*this); }

    template<typename ...args>
    DataType& element(const args&... index) { return (*this)(index...); }

    template<typename ...args>
    const DataType& element(const args&... index) const { return (*this)(index...); }

    template<typename ...args>
    DataType* element_ptr(const args&... index) { return &(*this)(index...); }

    template<typename ...args>
    const DataType* element_ptr(const args&... index) const { return &(*this)(index...); }

    RDM<rank,DataType>& operator+=(const RDM<rank,DataType>& o) { this->ax_plus_y(1.0, o); return *this; }
    RDM<rank,DataType>& operator-=(const RDM<rank,DataType>& o) { this->ax_plus_y(-1.0, o); return *this; }
    RDM<rank,DataType> operator+(const RDM<rank,DataType>& o) const { RDM<rank,DataType> out(*this); out += o; return out; }
    RDM<rank,DataType> operator-(const RDM<rank,DataType>& o) const { RDM<rank,DataType> out(*this); out -= o; return out; }

    void zero() { this->fill(0.0); }
    void allreduce() { mpi__->allreduce(data(), size()); }

    void ax_plus_y(const DataType a, const btas::TensorN<DataType,N>& o) { btas::axpy(a, o, *this); }
    void ax_plus_y(const DataType& a, const std::shared_ptr<const btas::TensorN<DataType,N>>& o) { this->ax_plus_y(a, *o); }

    void scale(const DataType& a) { btas::scal(a, *this); }
    size_t size() const { return this->storage().size(); }
    int norb() const { return this->extent(0); }

    // returns if this is natural orbitals - only for rank 1
    bool natural_orbitals() const {
      throw std::logic_error("RDM<N>::natural_orbitals() should not be called with N>1");
      return true;
    }

    std::shared_ptr<Matrix> rdm1_mat(const int nclosed, const bool all = true) const {
      throw std::logic_error("RDM<N>::rdm1_mat() should not be called with N>1");
      return nullptr;
    }

    std::shared_ptr<Matrix> rdm1_mat_tr(const int nclosed, const bool all = true) const {
      throw std::logic_error("RDM<N>::rdm1_mat_tr() should not be called with N>1");
      return nullptr;
    }

    void transform(std::shared_ptr<const Matrix> coeff) { throw std::logic_error("RDM<N>::transform() (N>3) not implemented yet"); }

    std::vector<DataType> diag() const {
      throw std::logic_error("RDM<N>::diag() should not be called with N>1");
      return std::vector<DataType>();
    }

    template<typename T = DataType,
             class = typename std::enable_if<std::is_same<T,std::complex<double>>::value>::type
            >
    std::shared_ptr<RDM<rank,double>> get_real_part() const {
      auto out = std::make_shared<RDM<rank,double>>(norb());
      auto i = out->begin();
      for (auto& d : *this)
        *i++ = std::real(d);
      return out;
    }

    template<typename T = DataType,
             class = typename std::enable_if<std::is_same<T,std::complex<double>>::value>::type
            >
    std::shared_ptr<RDM<rank,double>> get_imag_part() const {
      auto out = std::make_shared<RDM<rank,double>>(norb());
      auto i = out->begin();
      for (auto& d : *this)
        *i++ = std::imag(d);
      return out;
    }

    template<typename T = DataType,
             class = typename std::enable_if<std::is_same<T,std::complex<double>>::value>::type
            >
    std::shared_ptr<RDM<rank,std::complex<double>>> get_conjg() const {
      auto out = this->copy();
      for (auto& d : *out)
        d = std::conj(d);
      return out;
    }

    void print(const double thresh = 1.0e-3) const { throw std::logic_error("RDM<N>::print() (N>3) not implemented yet"); }
};

template <int rank>
using ZRDM = RDM<rank, std::complex<double>>;

template<> bool RDM<1,double>::natural_orbitals() const;
template<> std::vector<double> RDM<1,double>::diag() const;

template<> void RDM<1,double>::transform(std::shared_ptr<const Matrix> coeff);
template<> void RDM<2,double>::transform(std::shared_ptr<const Matrix> coeff);

template<> std::shared_ptr<Matrix> RDM<1,double>::rdm1_mat(const int nclosed, const bool all) const;
template<> std::shared_ptr<Matrix> RDM<1,double>::rdm1_mat_tr(const int nclosed, const bool all) const;

template<> void RDM<1,double>::print(const double thresh) const;
template<> void RDM<2,double>::print(const double thresh) const;
template<> void RDM<3,double>::print(const double thresh) const;
template<> void RDM<1,std::complex<double>>::print(const double thresh) const;
template<> void RDM<2,std::complex<double>>::print(const double thresh) const;
template<> void RDM<3,std::complex<double>>::print(const double thresh) const;

template<int rank>
using VecRDM = Vec<RDM<rank, double>>;

namespace detail {

template<int N, typename DataType, int M>
struct fill_in {
  fill_in() { }
  void operator()(std::shared_ptr<RDM<N,DataType>> out, std::shared_ptr<const RDM<N,DataType>> in, const std::array<size_t,2*N>& par, const size_t& norb,
                  const size_t offset1 = 0, const size_t offset2 = 0) {
    for (int i = 0; i != norb; ++i) {
      const size_t off1 = offset1*(2*norb) + par[M-1] + i;
      const size_t off2 = offset2*norb + i;
      fill_in<N,DataType,M-1>()(out, in, par, norb, off1, off2);
    }
  }
};

template<int N, typename DataType>
struct fill_in<N,DataType,2> {
  fill_in() { }
  void operator()(std::shared_ptr<RDM<N,DataType>> out, std::shared_ptr<const RDM<N,DataType>> in, const std::array<size_t,2*N>& par, const size_t& norb,
                  const size_t offset1 = 0, const size_t offset2 = 0) {
    for (int i = 0; i != norb; ++i) {
      const size_t off1 = offset1*(2*norb) + par[1] + i;
      const size_t off2 = offset2*norb + i;
      std::copy_n(in->data()+off2*norb, norb, out->data()+off1*2*norb+par[0]);
    }
  }
};
}


template<int N, typename DataType>
std::shared_ptr<RDM<N,DataType>> expand_kramers(std::shared_ptr<const Kramers<2*N,RDM<N,DataType>>> o, const size_t norb) {
//assert(!(o->begin()->second && o->begin()->second->norb() != norb));
  auto out = std::make_shared<RDM<N,DataType>>(2*norb);
  for (size_t i = 0; i != (1<<(2*N)); ++i) {
    auto data = o->get_data(i);
    if (data) {
      std::array<size_t,2*N> par;
      for (int j = 0; j != 2*N; ++j)
        par[2*N-1-j] = ((i>>j)&1)*norb;
      detail::fill_in<N,DataType,2*N>()(out, data, par, norb);
    }
  }
  return out;
}

}

#endif
