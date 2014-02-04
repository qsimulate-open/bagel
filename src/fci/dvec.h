//
// BAGEL - Parallel electron correlation program.
// Filename: dvec.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef BAGEL_FCI_DVEC_H
#define BAGEL_FCI_DVEC_H

#include <src/fci/civec.h>

// I forgot why I named this class "Dvector"...
// Basically Dvector is a vector of Civec, with some utility functions.
// It is used as intermediates in FCI, or eigenvectors in multistate CASSCF.

// I decided to make this class, in order to ensure that the date area
// in a FCI intermediate be consecutive so that I can use DGEMM etc.


// TODO The Dvector class is NOT yet flexible for Civectors with different Determinants objects.
// This can be easily done by modifing what follows.

namespace bagel {

template <typename DataType>
class Dvector {
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
    std::unique_ptr<DataType[]> data_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int version) const {
      // dvec_ is just an alias and therefore not serialized
      ar << det_ << lena_ << lenb_ << ij_ << boost::serialization::make_array(data(), size());
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int version) {
      ar >> det_ >> lena_ >> lenb_ >> ij_ >> dvec_;
      data_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      ar >> boost::serialization::make_array(data(), size());
      // make an alias and set to dvec_
      DataType* ptr = data();
      for (int i = 0; i != ij_; ++i, ptr += lena_*lenb_)
        dvec_.push_back(std::make_shared<Civector<DataType>>(det_, ptr));
    }

  public:
    Dvector() { }

    Dvector(std::shared_ptr<const Determinants> det, const size_t ij) : det_(det), lena_(det->lena()), lenb_(det->lenb()), ij_(ij) {
      // data should be in a contiguous area to call dgemm.
      data_ = std::unique_ptr<DataType[]>(new DataType[lenb_*lena_*ij_]);
      std::fill_n(data_.get(), lenb_*lena_*ij_, DataType(0.0));
      DataType* tmp = data_.get();
      for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_)
        dvec_.push_back(std::make_shared<Civector<DataType>>(det_, tmp));
    }

    Dvector(const Dvector<DataType>& o) : det_(o.det_), lena_(o.lena_), lenb_(o.lenb_), ij_(o.ij_) {
      data_ = std::unique_ptr<DataType[]>(new DataType[lena_*lenb_*ij_]);
      DataType* tmp = data_.get();
      for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_)
        dvec_.push_back(std::make_shared<Civector<DataType>>(det_, tmp));
      std::copy_n(o.data(), lena_*lenb_*ij_, data());
    }

    Dvector(std::shared_ptr<const Civector<DataType>> e, const size_t ij) : det_(e->det()), lena_(e->lena()), lenb_(e->lenb()), ij_(ij) {
      data_ = std::unique_ptr<DataType[]>(new DataType[lena_*lenb_*ij]);
      DataType* tmp = data();
      for (int i = 0; i != ij; ++i, tmp += lenb_*lena_) {
        auto c = std::make_shared<Civector<DataType>>(det_, tmp);
        std::copy_n(e->data(), lenb_*lena_, c->data());
        dvec_.push_back(c);
      }
    }

    // I think this is very confusiong... this is done this way in order not to delete Civec when Dvector is deleted.
    Dvector(std::shared_ptr<const Dvector<DataType>> o) : det_(o->det_), lena_(o->lena_), lenb_(o->lenb_), ij_(o->ij_) {
      for (int i = 0; i != ij_; ++i)
        dvec_.push_back(std::make_shared<Civector<DataType>>(*(o->data(i))));
    }

    Dvector(std::vector<std::shared_ptr<Civector<DataType>>> o) : det_(o.front()->det()), ij_(o.size()) {
      lena_ = det_->lena();
      lenb_ = det_->lenb();
      for (int i = 0; i != ij_; ++i)
        dvec_.push_back(std::make_shared<Civector<DataType>>(*(o.at(i))));
    }

    std::shared_ptr<const Determinants> det() const { return det_; }

    DataType* data() { return data_.get(); }
    const DataType* data() const { return data_.get(); }

    std::shared_ptr<Civector<DataType>>& data(const size_t i) { return dvec_[i]; }
    std::shared_ptr<const Civector<DataType>> data(const size_t i) const { return dvec_[i]; }
    void zero() { std::fill(data(), data()+lena_*lenb_*ij_, 0.0); }

    std::vector<std::shared_ptr<Civector<DataType>>>& dvec() { return dvec_; }
    const std::vector<std::shared_ptr<Civector<DataType>>>& dvec() const { return dvec_; }

    // returns a vector of Civec's which correspond to an unconverged state
    std::vector<std::shared_ptr<Civector<DataType>>> dvec(const std::vector<int>& conv) {
      std::vector<std::shared_ptr<Civector<DataType>>> out;
      auto c = conv.begin();
      for (auto& iter : dvec_) {
        if (*c++ == 0) out.push_back(iter);
        else out.push_back(std::shared_ptr<Civector<DataType>>());
      }
      return out;
    }
    std::vector<std::shared_ptr<const Civector<DataType>>> dvec(const std::vector<int>& conv) const {
      std::vector<std::shared_ptr<const Civector<DataType>>> out;
      auto c = conv.begin();
      for (auto& iter : dvec_) {
        if (*c++ == 0) out.push_back(iter);
        else out.push_back(std::shared_ptr<const Civector<DataType>>());
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

    void scale(const DataType& a) {
      std::for_each(dvec_.begin(), dvec_.end(), [&a](CiPtr p){ p->scale(a); });
    }

    Dvector& operator*=(const double& a) { scale(a); return *this; }

    std::shared_ptr<Dvector<DataType>> clone() const { return std::make_shared<Dvector<DataType>>(det_, ij_); }
    std::shared_ptr<Dvector<DataType>> copy() const { return std::make_shared<Dvector<DataType>>(*this); }

    // for double versions
    std::shared_ptr<Dvector<DataType>> spin() const { assert(false); return std::shared_ptr<Dvector<DataType>>(); }
    std::shared_ptr<Dvector<DataType>> spinflip(std::shared_ptr<const Determinants> det = std::shared_ptr<Determinants>()) const { assert(false); return std::shared_ptr<Dvector<DataType>>(); }
    std::shared_ptr<Dvector<DataType>> spin_lower(std::shared_ptr<const Determinants> det = std::shared_ptr<Determinants>()) const { assert(false); return std::shared_ptr<Dvector<DataType>>(); }
    std::shared_ptr<Dvector<DataType>> spin_raise(std::shared_ptr<const Determinants> det = std::shared_ptr<Determinants>()) const { assert(false); return std::shared_ptr<Dvector<DataType>>(); }

    std::shared_ptr<Dvector<DataType>> apply(const int orbital, const bool action, const int spin) const {
      std::vector<CiPtr> out;
      for (auto& i : dvec_) out.push_back(i->apply(orbital, action, spin));
      return std::make_shared<Dvector<DataType>>(out);
    }

    void orthog(std::shared_ptr<const Dvector<DataType>> o) {
      if (o->ij() != ij()) throw std::logic_error("Dvector<DataType>::orthog called inconsistently");
      std::transform(o->dvec_.begin(), o->dvec_.end(), dvec_.begin(), dvec_.begin(), [](CiPtr p, CiPtr q){ q->orthog(p); return q; });
    }

    void project_out(std::shared_ptr<const Dvector<DataType>> o) {
      if (o->ij() != ij()) throw std::logic_error("Dvec::project_out called inconsistently");
      auto j = o->dvec().begin();
      // simply project out each CI vector
      for (auto i = dvec().begin(); i != dvec().end(); ++i, ++j) (*i)->project_out(*j);
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

template<> std::shared_ptr<Dvector<double>> Dvector<double>::spin() const;
template<> std::shared_ptr<Dvector<double>> Dvector<double>::spinflip(std::shared_ptr<const Determinants> det) const;
template<> std::shared_ptr<Dvector<double>> Dvector<double>::spin_lower(std::shared_ptr<const Determinants> det) const;
template<> std::shared_ptr<Dvector<double>> Dvector<double>::spin_raise(std::shared_ptr<const Determinants> det) const;

using Dvec = Dvector<double>;
using ZDvec = Dvector<std::complex<double>>;

}

extern template class bagel::Dvector<double>;
extern template class bagel::Dvector<std::complex<double>>;

#endif
