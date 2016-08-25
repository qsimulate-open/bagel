//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dvec.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef BAGEL_FCI_DVEC_H
#define BAGEL_FCI_DVEC_H

#include <src/ci/fci/civec.h>

// I forgot why I named this class "Dvector"...
// Basically Dvector is a vector of Civec, with some utility functions.
// It is used as intermediates in FCI, or eigenvectors in multistate CASSCF.

// I decided to make this class, in order to ensure that the date area
// in a FCI intermediate be consecutive so that I can use DGEMM etc.


// TODO The Dvector class is NOT yet flexible for Civectors with different Determinants objects.
// This can be easily done by modifing what follows.

namespace bagel {

template <typename DataType>
class Dvector : public btas::Tensor3<DataType> {
  public:
    using data_type = DataType;
    using btas::Tensor3<DataType>::data;
    using btas::Tensor3<DataType>::begin;
    using btas::Tensor3<DataType>::cbegin;
    using btas::Tensor3<DataType>::end;
    using btas::Tensor3<DataType>::cend;

  // Useful in templates involving Dvectors
  public: using DetType = Determinants;
  public: using Ci    = Civector<DataType>;
  // only for use in lambdas
  private: using CiPtr = std::shared_ptr<Ci>;

  protected:
    // the determinant space where Dvector's are sitting
    mutable std::shared_ptr<const Determinants> det_;

    size_t lena_;
    size_t lenb_;
    // the size of the vector<shared_ptr<Civector<DataType>>>
    size_t ij_;
    std::vector<std::shared_ptr<Civector<DataType>>> dvec_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int version) const {
      ar << boost::serialization::base_object<btas::Tensor3<DataType>>(*this) << det_ << lena_ << lenb_ << ij_;
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int version) {
      ar >> boost::serialization::base_object<btas::Tensor3<DataType>>(*this) >> det_ >> lena_ >> lenb_ >> ij_;
      DataType* tmp = data();
      for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_)
        dvec_.push_back(std::make_shared<Civector<DataType>>(det_, tmp));
    }

  public:
    Dvector() { }

    Dvector(std::shared_ptr<const Determinants> det, const size_t ij) : btas::Tensor3<DataType>(det->lenb(), det->lena(), ij), det_(det), lena_(det->lena()), lenb_(det->lenb()), ij_(ij) {
      std::fill(begin(), end(), DataType(0.0));
      DataType* tmp = data();
      for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_)
        dvec_.push_back(std::make_shared<Civector<DataType>>(det_, tmp));
    }

    Dvector(const Dvector<DataType>& o) : btas::Tensor3<DataType>(o), det_(o.det_), lena_(o.lena_), lenb_(o.lenb_), ij_(o.ij_) {
      DataType* tmp = data();
      for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_)
        dvec_.push_back(std::make_shared<Civector<DataType>>(det_, tmp));
    }

    Dvector(const Dvector_base<Civector<DataType>>& o)
      : btas::Tensor3<DataType>(o.det()->lenb(), o.det()->lena(), o.ij()), det_(o.det()), lena_(o.det()->lena()), lenb_(o.det()->lenb()), ij_(o.ij()) {
      DataType* tmp = data();
      for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_)
        dvec_.push_back(std::make_shared<Civector<DataType>>(det_, tmp));
      auto iter = dvec_.begin();
      for (auto& inp : o.dvec())
        **iter++ = *inp;
    }

    Dvector(std::shared_ptr<const Civector<DataType>> e, const size_t ij) : btas::Tensor3<DataType>(e->lenb(), e->lena(), ij), det_(e->det()), lena_(e->lena()), lenb_(e->lenb()), ij_(ij) {
      DataType* tmp = data();
      for (int i = 0; i != ij; ++i, tmp += lenb_*lena_) {
        auto c = std::make_shared<Civector<DataType>>(det_, tmp);
        std::copy_n(e->data(), lenb_*lena_, c->data());
        dvec_.push_back(c);
      }
    }

    std::shared_ptr<const Determinants> det() const { return det_; }

    std::shared_ptr<Dvector> extract_state(const std::vector<int> input) const;

    std::shared_ptr<Civector<DataType>>& data(const size_t i) { return dvec_[i]; }
    std::shared_ptr<const Civector<DataType>> data(const size_t i) const { return dvec_[i]; }
    void zero() { std::fill(begin(), end(), 0.0); }

    std::vector<std::shared_ptr<Civector<DataType>>>& dvec() { return dvec_; }
    const std::vector<std::shared_ptr<Civector<DataType>>>& dvec() const { return dvec_; }

    // returns a vector of Civec's which correspond to an unconverged state
    std::vector<std::shared_ptr<Civector<DataType>>> dvec(const std::vector<int>& conv) {
      std::vector<std::shared_ptr<Civector<DataType>>> out;
      auto c = conv.begin();
      for (auto& i : dvec_) {
        if (*c++ == 0) out.push_back(i);
        else out.push_back(nullptr);
      }
      return out;
    }
    std::vector<std::shared_ptr<const Civector<DataType>>> dvec(const std::vector<int>& conv) const {
      std::vector<std::shared_ptr<const Civector<DataType>>> out;
      auto c = conv.begin();
      for (auto& i : dvec_) {
        if (*c++ == 0) out.push_back(i);
        else out.push_back(nullptr);
      }
      return out;
    }

    size_t lena() const { return lena_; }
    size_t lenb() const { return lenb_; }
    size_t ij() const { return ij_; }
    size_t size() const { return lena_*lenb_*ij_; }

    void set_det(std::shared_ptr<const Determinants> o) const {
      det_ = o;
      std::for_each(dvec_.begin(), dvec_.end(), [&o](CiPtr p){ p->set_det(o); });
    }

    // some functions for convenience
    DataType dot_product(const Dvector<DataType>& other) const {
      return std::inner_product(dvec_.begin(), dvec_.end(), other.dvec_.begin(), DataType(0.0), std::plus<DataType>(), [](CiPtr p, CiPtr q){ return p->dot_product(q); });
    }
    void ax_plus_y(const DataType a, std::shared_ptr<const Dvector<DataType>> other) {
      ax_plus_y(a, *other);
    }
    void ax_plus_y(const DataType a, const Dvector<DataType>& other) {
      std::transform(other.dvec_.begin(), other.dvec_.end(), dvec_.begin(), dvec_.begin(), [&a](CiPtr p, CiPtr q){ q->ax_plus_y(a, p); return q; });
    }
    Dvector<DataType>& operator+=(const Dvector<DataType>& o) { ax_plus_y(1.0, o); return *this; }
    Dvector<DataType>& operator-=(const Dvector<DataType>& o) { ax_plus_y(-1.0, o); return *this; }

    Dvector<DataType> operator+(const Dvector<DataType>& o) const { Dvector<DataType> out(*this); return out += o; }
    Dvector<DataType> operator-(const Dvector<DataType>& o) const { Dvector<DataType> out(*this); return out -= o; }

    Dvector<DataType>& operator/=(const Dvector<DataType>& o) {
      assert(dvec().size() == o.dvec().size());
      std::transform(o.dvec_.begin(), o.dvec_.end(), dvec_.begin(), dvec_.begin(), [](CiPtr p, CiPtr q){ *q / *p; return q; });
      return *this;
    }
    Dvector<DataType> operator/(const Dvector<DataType>& o) const {
      Dvector<DataType> out(*this);
      out /= o;
      return out;
    }

    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    double rms() const { return norm() / std::sqrt(size()); }

    void scale(const DataType& a) {
      std::for_each(dvec_.begin(), dvec_.end(), [&a](CiPtr p){ p->scale(a); });
    }

    Dvector& operator*=(const double& a) { scale(a); return *this; }

    std::shared_ptr<Dvector<DataType>> clone() const { return std::make_shared<Dvector<DataType>>(det_, ij_); }
    std::shared_ptr<Dvector<DataType>> copy() const { return std::make_shared<Dvector<DataType>>(*this); }

    void orthog(std::shared_ptr<const Dvector<DataType>> o) {
      if (o->ij() != ij()) throw std::logic_error("Dvector<DataType>::orthog called inconsistently");
      std::transform(o->dvec_.begin(), o->dvec_.end(), dvec_.begin(), dvec_.begin(), [](CiPtr p, CiPtr q){ q->orthog(p); return q; });
    }

    void project_out(std::shared_ptr<const Dvector<DataType>> o) {
      for (auto& i : dvec())
        for (auto& j : o->dvec())
          i->project_out(j);
    }

    void synchronize() {
      for (auto& i : dvec_)
        i->synchronize();
    }

    void print(const double thresh = 0.05) const {
      int j = 0;
      for (auto& iter : dvec_) {
        std::cout << std::endl << "     * ci vector, state " << std::setw(3) << j++;
        if (typeid(DataType) == typeid(double)) std::cout << ", <S^2> = " << std::setw(6) << std::setprecision(4) << iter->spin_expectation();
        std::cout << std::endl;
        iter->print(thresh);
      }
    }
};

using Dvec = Dvector<double>;
using ZDvec = Dvector<std::complex<double>>;

}

extern template class bagel::Dvector<double>;
extern template class bagel::Dvector<std::complex<double>>;

#endif
