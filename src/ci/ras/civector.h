//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/civector.h
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


#ifndef __BAGEL_RAS_RASCIVECTOR_H
#define __BAGEL_RAS_RASCIVECTOR_H

#include <src/ci/ras/civector_impl.h>

namespace bagel {

template <typename DataType> class RASCivecView_;

template <typename DataType>
class RASCivector : public RASCivector_impl<DataType> {
  public:
    using DetType = RASDeterminants;
    using RBlock = RASBlock<DataType>;
    using LocalizedType = std::true_type;
    using RASCivector_base<RBlock>::det;
    using RASCivector_base<RBlock>::size;

  protected:
    using RASCivector_base<RASBlock<DataType>>::blocks_;

    std::unique_ptr<DataType[]> data_;

    std::shared_ptr<RASCivector_impl<DataType>> spin_() const { return spin(); }

  public:
    RASCivector(std::shared_ptr<const RASDeterminants> det);

    RASCivector(const RASCivector<DataType>& o) : RASCivector(o.det()) { std::copy_n(o.data(), size(), data_.get()); }
    RASCivector(RASCivector<DataType>&& o) : RASCivector_impl<DataType>(o.det()) { blocks_ = std::move(o.blocks_); }

    RASCivector(const RASCivecView_<DataType>& o) : RASCivector(o.det()) { std::copy_n(o.data(), size(), data_.get()); }

    RASCivector(std::shared_ptr<const RASCivector<DataType>> o) : RASCivector(o->det()) { std::copy_n(o->data(), size(), data_.get()); }
    RASCivector(std::shared_ptr<const RASCivecView_<DataType>>& o) : RASCivector(o->det()) { std::copy_n(o->data(), size(), data_.get()); }

    // Move assignment
    RASCivector<DataType>& operator=(RASCivector<DataType>&& o);

    DataType* data() override { return data_.get(); }
    const DataType* data() const override { return data_.get(); }

    std::shared_ptr<RASCivector<DataType>> clone() const { return std::make_shared<RASCivector<DataType>>(det()); }
    std::shared_ptr<RASCivector<DataType>> copy() const  { return std::make_shared<RASCivector<DataType>>(*this); }

    // Spin functions are only implememted as specialized functions for double (see civec.cc)
    std::shared_ptr<RASCivector<DataType>> spin() const;
    std::shared_ptr<RASCivector<DataType>> spin_lower(std::shared_ptr<const RASDeterminants> target_det = nullptr) const;
    std::shared_ptr<RASCivector<DataType>> spin_raise(std::shared_ptr<const RASDeterminants> target_det = nullptr) const;

    std::shared_ptr<RASCivector<DataType>> transpose(std::shared_ptr<const RASDeterminants> det = nullptr) const {
      std::shared_ptr<RASCivector<DataType>> out;
      this->transpose_impl(out, det);
      return out;
    }

    std::shared_ptr<RASCivector<DataType>> apply(const int orbital, const bool action, const bool spin) const;
};

template<> std::shared_ptr<RASCivector<double>> RASCivector<double>::spin() const; // returns S^2 | civec >
template<> std::shared_ptr<RASCivector<double>> RASCivector<double>::spin_lower(std::shared_ptr<const RASDeterminants>) const; // S_-
template<> std::shared_ptr<RASCivector<double>> RASCivector<double>::spin_raise(std::shared_ptr<const RASDeterminants>) const; // S_+

using RASCivec = RASCivector<double>;
using RASDvec  = Dvector_base<RASCivec>;

/*********************************************************************************************************************************/

template <typename DataType>
class RASCivecView_ : public RASCivector_impl<DataType> {
  public:
    using DetType = RASDeterminants;
    using RBlock = RASBlock<DataType>;
    using LocalizedType = std::true_type;
    using RASCivector_base<RBlock>::det;
    using RASCivector_base<RBlock>::size;

  protected:
    using RASCivector_base<RASBlock<DataType>>::blocks_;

    double* const data_ptr_;
    bool can_write_;

    std::shared_ptr<RASCivector_impl<DataType>> spin_() const { return spin(); }

  public:
    RASCivecView_(std::shared_ptr<const RASDeterminants> det, double* const data);
    RASCivecView_(std::shared_ptr<const RASDeterminants> det, const double* const data)
     : RASCivecView_(det, const_cast<DataType*>(data)) { can_write_ = false; }

    DataType* data() override { assert(can_write_); return data_ptr_; }
    const DataType* data() const override { return data_ptr_; }

    RASCivecView_(RASCivector<DataType>& o) : RASCivecView_(o.det(), o.data()) {}
    RASCivecView_(const RASCivector<DataType>& o) : RASCivecView_(o.det(), o.data()) {}
    RASCivecView_(const RASCivecView_<DataType>& o) : RASCivecView_(o.det(), o.data()) { can_write_ = o.can_write_; }

    // Spin functions are only implememted as specialized functions for double (see civec.cc)
    std::shared_ptr<RASCivector<DataType>> spin() const;
    std::shared_ptr<RASCivector<DataType>> spin_lower(std::shared_ptr<const RASDeterminants> target_det = nullptr) const;
    std::shared_ptr<RASCivector<DataType>> spin_raise(std::shared_ptr<const RASDeterminants> target_det = nullptr) const;

    std::shared_ptr<RASCivector<DataType>> transpose(std::shared_ptr<const RASDeterminants> det = nullptr) const {
      std::shared_ptr<RASCivector<DataType>> out;
      this->transpose_impl(out, det);
      return out;
    }
};

template<> std::shared_ptr<bagel::RASCivector<double>> RASCivecView_<double>::spin() const; // returns S^2 | civec >
template<> std::shared_ptr<RASCivector<double>> RASCivecView_<double>::spin_lower(std::shared_ptr<const RASDeterminants>) const; // S_-
template<> std::shared_ptr<RASCivector<double>> RASCivecView_<double>::spin_raise(std::shared_ptr<const RASDeterminants>) const; // S_+

using RASCivecView = RASCivecView_<double>;

extern template class RASCivector<double>;
extern template class RASCivecView_<double>;

}

#endif
