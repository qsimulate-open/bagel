//
// BAGEL - Parallel electron correlation program.
// Filename: fci/dvector.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef BAGEL_RAS_DVECTOR_H
#define BAGEL_RAS_DVECTOR_H

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
  public: using DetType = typename CiType::DetType;
  public: using Ci = CiType;
  // only for use in lambdas
  using CiPtr = std::shared_ptr<CiType>;

  protected:
    // the determinant space where Dvector_base's are sitting
    mutable std::shared_ptr<const DetType> det_;

    // the size of the vector<CiPtr>
    size_t ij_;

    std::vector<CiPtr> dvec_;

  public:
    Dvector_base(std::shared_ptr<const DetType> det, const size_t ij) : det_(det), ij_(ij) {
      dvec_.clear();
      for (int i = 0; i < ij_; ++i) dvec_.push_back( std::make_shared<CiType>(det_) );
    }

    Dvector_base(const Dvector_base<CiType>& o) : det_(o.det_), ij_(o.ij_) {
      dvec_.clear();
      for ( auto& ivec : o.dvec() ) dvec_.push_back( std::make_shared<CiType>(*ivec) );
    }
    Dvector_base(std::shared_ptr<const Dvector_base<CiType>> o) : Dvector_base<CiType>(*o) {}

    Dvector_base(std::vector<CiPtr> o) : det_(o.front()->det()), ij_(o.size()) {
      for (auto& ivec : o) dvec_.push_back( std::make_shared<CiType>(*ivec) );
    }

    std::shared_ptr<const DetType> det() const { return det_; }

    CiPtr& data(const size_t i) { return dvec_[i]; }
    CiPtr data(const size_t i) const { return dvec_[i]; }

    std::vector<CiPtr>& dvec() { return dvec_; }
    std::vector<CiPtr> dvec() const { return dvec_; }

    size_t ij() const { return ij_; }
    size_t size() const { return ij_ * det_->size(); }

    void set_det(std::shared_ptr<const DetType> o) const {
      det_ = o;
      std::for_each(dvec_.begin(), dvec_.end(), [&o](CiPtr p) { p->set_det(o); });
    }

    std::shared_ptr<Dvector_base<CiType>> clone() const { return std::make_shared<Dvector_base<CiType>>(det_, ij_); }
    std::shared_ptr<Dvector_base<CiType>> copy() const { return std::make_shared<Dvector_base<CiType>>(*this); }

    // for MEH
    std::shared_ptr<Dvector_base<CiType>> apply(const int orbital, const bool action, const int spin) const {
      std::vector<CiPtr> out;
      for (auto& i : dvec_) out.push_back(i->apply(orbital, action, spin));
      return std::make_shared<Dvector_base<CiType>>(out);
    }


    // will fail for non-double DataTypes
    std::shared_ptr<Dvector_base<CiType>> spin() const {
      std::vector<CiPtr> out;
      for (auto& i : dvec_) out.push_back( i->spin() );
      return std::make_shared<Dvector_base<CiType>>(out);
    }

    std::shared_ptr<Dvector_base<CiType>> transpose(std::shared_ptr<const DetType> det = std::shared_ptr<DetType>()) const {
      if (!det) det = det_->transpose();
      std::vector<CiPtr> out;
      for (auto& i : dvec_) out.push_back( i->transpose(det) );
      return std::make_shared<Dvector_base<CiType>>(out);
    }

    std::shared_ptr<Dvector_base<CiType>> spin_lower(std::shared_ptr<const DetType> det = std::shared_ptr<DetType>()) const {
      if (!det) det = det_->clone(det_->nelea() - 1, det_->neleb() + 1);
      std::vector<CiPtr> out;
      for (auto& i : dvec_) out.push_back( i->spin_lower(det) );
      return std::make_shared<Dvector_base<CiType>>(out);
    }

    std::shared_ptr<Dvector_base<CiType>> spin_raise(std::shared_ptr<const DetType> det = std::shared_ptr<DetType>()) const {
      if (!det) det = det_->clone(det_->nelea() + 1, det_->neleb() - 1);
      std::vector<CiPtr> out;
      for (auto& i : dvec_) out.push_back( i->spin_raise(det) );
      return std::make_shared<Dvector_base<CiType>>(out);
    }

    void orthog(std::shared_ptr<const Dvector_base<CiType>> o) {
      if (o->ij() != ij()) throw std::logic_error("Dvector_base<CiType>::orthog called inconsistently");
      std::transform(o->dvec_.begin(), o->dvec_.end(), dvec_.begin(), dvec_.begin(), [](CiPtr p, CiPtr q){ q->orthog(p); return q; });
    }

    void project_out(std::shared_ptr<const Dvector_base<CiType>> o) {
      if (o->ij() != ij()) throw std::logic_error("Dvector_base<CiType>::project_out called inconsistently");
      auto j = o->dvec().begin();
      // simply project out each CI vector
      for (auto i = dvec().begin(); i != dvec().end(); ++i, ++j) (*i)->project_out(*j);
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
};

}

#endif
