//
// BAGEL - Parallel electron correlation program.
// Filename: ras/dvec.h
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


#ifndef BAGEL_RAS_DVEC_H
#define BAGEL_RAS_DVEC_H

#include <src/ras/civec.h>

// This class is purely for convenience and compliance with FCI Dvec.
// Unlike the FCI version though, this is basically just a wrapper around a vector

namespace bagel {

template <typename DataType>
class RASDvector {
  // only for use in lambdas
  using Ci    = RASCivector<DataType>;
  using CiPtr = std::shared_ptr<Ci>;

  protected:
    // the determinant space where Dvector's are sitting
    mutable std::shared_ptr<const RASDeterminants> det_;

    // the size of the vector<CiPtr>
    size_t ij_;

    std::vector<CiPtr> dvec_;

  public:
    RASDvector(std::shared_ptr<const RASDeterminants> det, const size_t ij) : det_(det), ij_(ij) {
      dvec_.clear();
      for (int i = 0; i < ij_; ++i) dvec_.push_back( std::make_shared<Ci>(det_) );
    }

    RASDvector(const RASDvector<DataType>& o) : det_(o.det_), ij_(o.ij_) {
      dvec_.clear();
      for ( auto& ivec : o.dvec() ) dvec_.push_back( std::make_shared<Ci>(*ivec) );
    }

    RASDvector(std::vector<CiPtr> o) : det_(o.front()->det()), ij_(o.size()) {
      for (int i = 0; i != ij_; ++i)
        dvec_.push_back(std::make_shared<Civector<DataType>>(*(o.at(i))));
    }

    std::shared_ptr<const RASDeterminants> det() const { return det_; }

    CiPtr& data(const size_t i) { return dvec_[i]; }
    CiPtr data(const size_t i) const { return dvec_[i]; }

    std::vector<CiPtr>& dvec() { return dvec_; }
    std::vector<CiPtr> dvec() const { return dvec_; }

    size_t ij() const { return ij_; }
    size_t size() const { return ij_ * det->size(); }

    void set_det(std::shared_ptr<const RASDeterminants> o) const {
      det_ = o;
      std::for_each(dvec_.begin(), dvec_.end(), [&o](CiPtr p){ p->set_det(o); });
    }

    std::shared_ptr<RASDvector<DataType>> clone() const { return std::make_shared<RASDvector<DataType>>(det_, ij_); }
    std::shared_ptr<RASDvector<DataType>> copy() const { return std::make_shared<RASDvector<DataType>>(*this); }

    // will fail for non-double DataTypes
    std::shared_ptr<RASDvector<DataType>> spin() const {
      std::vector<CiPtr> out;
      for (auto& i : dvec_) out.push_back( i->spin() );
      return std::make_shared<RASDvector<DataType>>(out);
    }

    std::shared_ptr<RASDvector<DataType>> transpose(std::shared_ptr<const RASDeterminants> det = std::shared_ptr<RASDeterminants>()) const {
      if (!det) det = det_->transpose();
      std::vector<CiPtr> out;
      for (auto& i : dvec_) out.push_back( i->transpose() );
      return std::make_shared<RASDvector<DataType>>(out);
    }

    std::shared_ptr<RASDvector<DataType>> spin_lower(std::shared_ptr<const RASDeterminants> det = std::shared_ptr<RASDeterminants>()) const {
      if (!det) det = det_->clone(det_->nelea() - 1, det_->neleb() + 1);
      std::vector<CiPtr> out;
      for (auto& i : dvec_) out.push_back( i->spin_lower );
      return std::make_shared<RASDvector<DataType>>(out);
    }

    std::shared_ptr<RASDvector<DataType>> spin_raise(std::shared_ptr<const RASDeterminants> det = std::shared_ptr<RASDeterminants>()) const {
      if (!det) det = det_->clone(det_->nelea() + 1, det_->neleb() - 1);
      std::vector<CiPtr> out;
      for (auto& i : dvec_) out.push_back( i->spin_raise );
      return std::make_shared<RASDvector<DataType>>(out);
    }

    void orthog(std::shared_ptr<const RASDvector<DataType>> o) {
      if (o->ij() != ij()) throw std::logic_error("RASDvector<DataType>::orthog called inconsistently");
      std::transform(o->dvec_.begin(), o->dvec_.end(), dvec_.begin(), dvec_.begin(), [](CiPtr p, CiPtr q){ q->orthog(p); return q; });
    }

    void project_out(std::shared_ptr<const RASDvector<DataType>> o) {
      if (o->ij() != ij()) throw std::logic_error("RASDvec::project_out called inconsistently");
      auto j = o->dvec().begin();
      // simply project out each CI vector
      for (auto i = dvec().begin(); i != dvec().end(); ++i, ++j) (*i)->project_out(*j);
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

#endif
