//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci/dvector.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef BAGEL_FCI_DVECTOR_BASE_H
#define BAGEL_FCI_DVECTOR_BASE_H

#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

// This version of Dvector_base is mainly just a wrapper around some vectors

namespace bagel {

template <typename CiType>
class Dvector_base {
  // used for template magicking
  public:
    using DetType = typename CiType::DetType;
    using Ci = CiType;
    using CiPtr = std::shared_ptr<CiType>;
    using CiConstPtr = std::shared_ptr<const CiType>;

  protected:
    // the determinant space where Dvector_base's are sitting
    mutable std::shared_ptr<const DetType> det_;

    // the size of the vector<CiPtr>
    size_t ij_;

    std::vector<CiPtr> dvec_;

  public:
    Dvector_base(std::shared_ptr<const DetType> det, const size_t ij) : det_(det), ij_(ij) {
      dvec_.clear();
      for (int i = 0; i < ij_; ++i) dvec_.push_back(std::make_shared<CiType>(det_));
    }

    Dvector_base(const Dvector_base<CiType>& o) : det_(o.det_), ij_(o.ij_) {
      for (auto& ivec : o.dvec()) dvec_.push_back(std::make_shared<CiType>(*ivec));
    }

    Dvector_base(const std::vector<CiPtr> o) : det_(o.front()->det()), ij_(o.size()) {
      for (auto& ivec : o) dvec_.push_back(ivec->copy());
    }

    template <class T>
    Dvector_base(std::shared_ptr<const Dvector_base<T>> o) : det_(o->det()), ij_(o->ij()) {
      for (int i = 0; i < ij_; ++i) dvec_.push_back(std::make_shared<CiType>(o->data(i)));
    }

    void zero() { std::for_each(dvec_.begin(), dvec_.end(), [](CiPtr p) { p->zero(); }); }

    template <typename DataType>
    void ax_plus_y(const DataType& a, const Dvector_base<CiType>& o) {
      assert(o.ij() == ij());
      auto oiter = o.dvec().begin();
      for (auto& i : dvec_)
        i->ax_plus_y(a, *oiter++);
    }

    std::shared_ptr<const DetType> det() const { return det_; }

    CiPtr& data(const size_t i) { return dvec_[i]; }
    CiConstPtr data(const size_t i) const { return dvec_[i]; }

    std::vector<CiPtr>& dvec() { return dvec_; }
    const std::vector<CiPtr>& dvec() const { return dvec_; }

    // returns a vector of CiType's which correspond to an unconverged state
    std::vector<CiConstPtr> dvec(const std::vector<int>& conv) const {
      std::vector<CiConstPtr> out;
      auto c = conv.begin();
      for (auto& i : dvec_) {
        if (*c++ == 0) out.push_back(i);
        else out.push_back(nullptr);
      }
      return out;
    }

    size_t ij() const { return ij_; }
    size_t size() const { return ij_ * det_->size(); }

    void set_det(std::shared_ptr<const DetType> o) const {
      det_ = o;
      std::for_each(dvec_.begin(), dvec_.end(), [&o](CiPtr p) { p->set_det(o); });
    }

    std::shared_ptr<Dvector_base<CiType>> clone() const { return std::make_shared<Dvector_base<CiType>>(det_, ij_); }
    std::shared_ptr<Dvector_base<CiType>> copy() const { return std::make_shared<Dvector_base<CiType>>(*this); }

    // for ASD
    std::shared_ptr<Dvector_base<CiType>> apply(const int orbital, const bool action, const bool spin) const {
      std::vector<std::shared_ptr<CiType>> out;
      for (auto& i : dvec_) out.push_back( i->apply(orbital, action, spin) );
      return std::make_shared<Dvector_base<CiType>>(out);
    }

    std::shared_ptr<Dvector_base<CiType>> apply_and_allocate(const bool action, const bool spin) const {
      std::shared_ptr<const DetType> source_det = this->det();
      std::shared_ptr<const DetType> target_det = spin ? (action ? source_det->addalpha() : source_det->remalpha()) :
                                                         (action ? source_det->addbeta() : source_det->rembeta());
      return std::make_shared<Dvector_base<CiType>>(target_det, this->ij());
    }

    void apply_and_fill(std::shared_ptr<const Dvector_base<CiType>> source_dvec, const int orbital, const bool action, const bool spin) { assert(false); }

    // will fail for non-double DataTypes
    std::shared_ptr<Dvector_base<CiType>> spin() const {
      return form_from_each([] (std::shared_ptr<const CiType> cc) { return cc->spin(); }, det(), typename CiType::LocalizedType());
    }

    std::shared_ptr<Dvector_base<CiType>> transpose(std::shared_ptr<const DetType> det = nullptr) const {
      if (!det) det = det_->transpose();
      return transpose_impl(det, typename CiType::LocalizedType());
    }

    std::shared_ptr<Dvector_base<CiType>> spin_lower(std::shared_ptr<const DetType> det = nullptr) const {
      if (!det) det = det_->clone(det_->nelea() - 1, det_->neleb() + 1);
      return form_from_each([det] (std::shared_ptr<const CiType> cc) { return cc->spin_lower(det); }, det, typename CiType::LocalizedType());
    }

    std::shared_ptr<Dvector_base<CiType>> spin_raise(std::shared_ptr<const DetType> det = nullptr) const {
      if (!det) det = det_->clone(det_->nelea() + 1, det_->neleb() - 1);
      return form_from_each([det] (std::shared_ptr<const CiType> cc) { return cc->spin_raise(det); }, det, typename CiType::LocalizedType());
    }

    void orthog(std::shared_ptr<const Dvector_base<CiType>> o) {
      if (o->ij() != ij()) throw std::logic_error("Dvector_base<CiType>::orthog called inconsistently");
      std::transform(o->dvec_.begin(), o->dvec_.end(), dvec_.begin(), dvec_.begin(), [](CiPtr p, CiPtr q){ q->orthog(p); return q; });
    }

    void project_out(std::shared_ptr<const Dvector_base<CiType>> o) {
      if (o->ij() != ij()) throw std::logic_error("Dvector_base<CiType>::project_out called inconsistently");
      auto j = o->dvec().begin();
      for (auto& i : dvec())
        i->project_out(*j++);
    }

    void project_out_all(std::shared_ptr<const Dvector_base<CiType>> o) {
      for (auto& i : dvec())
        for (auto& j : o->dvec())
          i->project_out(j);
    }

    void print(const double thresh = 0.05) const {
      int j = 0;
      for (auto& iter : dvec_) {
        std::cout << std::endl << "     * ci vector, state " << std::setw(3) << j++;
        std::cout << ", <S^2> = " << std::setw(6) << std::setprecision(4) << iter->spin_expectation();
        std::cout << std::endl;
        iter->print(thresh);
      }
    }

  private:
    // applies given function to each civec and returns a Dvector built from the function's return values
    // supplied determinant is used to construct blank Civecs
    // function must be of type
    //   std::shared_ptr<CiType> function(std::shared_ptr<const CiType>);
    // local version
    template <class Func>
    std::shared_ptr<Dvector_base<CiType>> form_from_each(Func func, std::shared_ptr<const DetType> d, /*local*/ std::true_type) const {
      std::vector<CiPtr> out;
      for (int i = 0; i < ij_; ++i) {
#ifdef HAVE_MPI_H
        if ((i % mpi__->size()) == mpi__->rank())
          out.push_back(func(data(i)));
        else
          out.push_back(std::make_shared<CiType>(d));
#else
        out.push_back(func(data(i)));
#endif
      }
#ifdef HAVE_MPI_H
      for (int i = 0; i < ij_; ++i) {
        out[i]->synchronize(i % mpi__->size());
      }
#endif
      return std::make_shared<Dvector_base<CiType>>(out);
    }

    // dist version
    template <class Func>
    std::shared_ptr<Dvector_base<CiType>> form_from_each(Func func, std::shared_ptr<const DetType> d, /*local*/ std::false_type) const {
      std::vector<CiPtr> out;
      for (int i = 0; i < ij_; ++i)
        out.push_back( func(data(i)) );
      return std::make_shared<Dvector_base<CiType>>(out);
    }

    // local transpose
    std::shared_ptr<Dvector_base<CiType>> transpose_impl(std::shared_ptr<const DetType> det, std::true_type) const {
      return form_from_each( [&det] (std::shared_ptr<const CiType> c) { return c->transpose(det); }, det, std::true_type() );
    }

    // dist transpose
    std::shared_ptr<Dvector_base<CiType>> transpose_impl(std::shared_ptr<const DetType> det, std::false_type) const {
      std::vector<CiPtr> out;
      for (auto& c : dvec_) {
        CiPtr tmp = c->transpose(det);
        tmp->transpose_wait();
        out.push_back(tmp);
      }
      return std::make_shared<CiType>(out);
    }

};

}

#endif
