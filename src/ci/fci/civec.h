//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: civec.h
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


#ifndef BAGEL_FCI_CIVEC_H
#define BAGEL_FCI_CIVEC_H

#include <list>
#include <src/util/math/algo.h>
#include <src/util/f77.h>
#include <src/util/parallel/staticdist.h>
#include <src/ci/fci/determinants.h>
#include <src/ci/fci/dvector_base.h>

namespace bagel {

template<typename DataType>
class Civector {
  public: using DetType = Determinants; // used to automatically determine type for Determinants object in templates
  public: using LocalizedType = std::true_type;

  protected:
    // The determinant space in which this Civec object is defined
    mutable std::shared_ptr<const Determinants> det_;

    size_t lena_;
    size_t lenb_;

    // !!CAUTION!!
    // cc is formated so that B runs first.
    // Also, cc_ can be null if this is constructed by Dvec.
    std::unique_ptr<DataType[]> cc_;

    DataType* cc_ptr_;

    DataType& cc(const size_t& i) { return *(cc_ptr_+i); }
    const DataType& cc(const size_t& i) const { return *(cc_ptr_+i); }
    DataType* cc() { return cc_ptr_; }
    const DataType* cc() const { return cc_ptr_; }

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      if (!cc_.get())
        throw std::logic_error("illegal call of Civector<T>::save");
      ar << det_ << lena_ << lenb_ << make_array(cc(), size());
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> det_ >> lena_ >> lenb_;
      cc_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      cc_ptr_ = cc_.get();
      ar >> make_array(cc(), size());
    }

  public:
    Civector() { }
    Civector(std::shared_ptr<const Determinants> det);
    // constructor that is called by Dvec.
    Civector(std::shared_ptr<const Determinants> det, DataType* din_);
    // copy constructor
    Civector(const Civector<DataType>& o);
    // this is not a copy constructor.
    Civector(std::shared_ptr<Civector<DataType>> o, std::shared_ptr<const Determinants> det);

    std::shared_ptr<Civector<DataType>> clone() const { return std::make_shared<Civector<DataType>>(det_); }
    std::shared_ptr<Civector<DataType>> copy() const { return std::make_shared<Civector<DataType>>(*this); }

    DataType* data() { return cc(); }
    DataType& element(size_t i, size_t j) { return cc(i+j*lenb_); }
    DataType* element_ptr(size_t i, size_t j) { return cc()+i+j*lenb_; }

    DataType& data(const size_t& i) { return cc(i); }
    const DataType& data(const size_t& i) const { return cc(i); }

    const DataType* data() const { return cc(); }
    const DataType* element_ptr(size_t i, size_t j) const { return cc()+i+j*lenb_; }

    std::shared_ptr<const Determinants> det() const { return det_; }
    void set_det(std::shared_ptr<const Determinants> o) const { det_ = o; }

    void zero() { std::fill(cc(), cc()+lena_*lenb_, DataType(0.0)); }

    size_t size() const { return lena_*lenb_; }

    std::shared_ptr<Civector<DataType>> transpose(std::shared_ptr<const Determinants> det = nullptr) const;

    size_t lena() const { return lena_; }
    size_t lenb() const { return lenb_; }
    size_t asize() const { return lena_; }

    // some functions for convenience
    void ax_plus_y(DataType a, const Civector<DataType>& other) {
      assert((lena_ == other.lena_) && (lenb_ == other.lenb_));
      blas::ax_plus_y_n(a, other.data(), size(), cc());
    }
    void ax_plus_y(DataType a, std::shared_ptr<const Civector> other) { ax_plus_y(a, *other); }

    DataType dot_product(const Civector<DataType>& other) const {
      assert((lena_ == other.lena_) && (lenb_ == other.lenb_));
      return blas::dot_product(cc(), size(), other.data());
    }
    DataType dot_product(std::shared_ptr<const Civector> other) const { return dot_product(*other); }

    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    double variance() const { return detail::real(dot_product(*this)) / size(); }
    double rms() const { return std::sqrt(variance()); }

    void scale(const DataType a) { blas::scale_n(a, cc(), size()); }

    // Spin functions are only implememted as specialized functions for double (see civec.cc)
    DataType spin_expectation() const { // returns < S^2 >
      std::shared_ptr<Civector<DataType>> S2 = spin();
      return dot_product(*S2);
    }
    std::shared_ptr<Civector<DataType>> spin() const;

    // S_- = \sum_i i_beta^\dagger i_alpha
    std::shared_ptr<Civector<DataType>> spin_lower(std::shared_ptr<const Determinants> target_det = nullptr) const;

    // S_+ = \sum_i i_alpha^\dagger i_beta
    std::shared_ptr<Civector<DataType>> spin_raise(std::shared_ptr<const Determinants> target_det = nullptr) const;

    void spin_decontaminate(const double thresh = 1.0e-12);

    std::shared_ptr<Civector<DataType>> apply(const int orbital, const bool action, const bool spin) const;

    Civector<DataType>& operator*=(const double& a) { scale(a); return *this; }
    Civector<DataType>& operator+=(const double& a) { std::for_each(cc(), cc()+size(), [&a](DataType& p){ p += a; }); return *this; }
    Civector<DataType>& operator-=(const double& a) { std::for_each(cc(), cc()+size(), [&a](DataType& p){ p -= a; }); return *this; }

    Civector<DataType>& operator=(const Civector<DataType>& o) {
      assert(det()->lena() == o.det()->lena() && det()->lenb() == o.det()->lenb());
      std::copy_n(o.cc(), size(), cc());
      return *this;
    }
    Civector<DataType>& operator+=(const Civector<DataType>& o) { ax_plus_y( 1.0, o); return *this; }
    Civector<DataType>& operator-=(const Civector<DataType>& o) { ax_plus_y(-1.0, o); return *this; }
    Civector<DataType>& operator/=(const Civector<DataType>& o) {
      for (size_t i = 0; i != size(); ++i)
        data(i) /= o.data(i);
      return *this;
    }
    Civector<DataType> operator/(const Civector<DataType>& o) const {
      Civector<DataType> out(*this);
      out /= o;
      return out;
    }

    // assumes that Civec's in c are already orthogonal with each other.
    // returns scaling factor (see implementation)

    double orthog(std::list<std::shared_ptr<const Civector<DataType>>> c) {
      for (auto& iter : c)
        project_out(iter);
      return normalize();
    }

    double orthog(std::shared_ptr<const Civector<DataType>> o) {
      return orthog(std::list<std::shared_ptr<const Civector<DataType>>>{o});
    }

    double normalize() {
      const double norm = this->norm();
      const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
      scale(scal);
      return 1.0/scal;
    }

    void project_out(std::shared_ptr<const Civector<DataType>> o) { ax_plus_y(-detail::conj(dot_product(*o)), *o); }

    void print(const double thr, const bool sort = true) const;

    void synchronize(const int root = 0) {
      mpi__->broadcast(cc_ptr_, size(), root);
    }
};

template<> void Civector<double>::spin_decontaminate(const double thresh);
template<> void Civector<std::complex<double>>::spin_decontaminate(const double thresh);

using Civec = Civector<double>;
using ZCivec = Civector<std::complex<double>>;

template <>
void Dvector_base<Civec>::apply_and_fill(std::shared_ptr<const Dvector_base<Civec>> source_dvec, const int orbital, const bool action, const bool spin);

using CASDvec = Dvector_base<Civec>;

}

extern template class bagel::Civector<double>;
extern template class bagel::Civector<std::complex<double>>;

#endif
