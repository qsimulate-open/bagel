//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: tensor.h
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


#ifndef __SRC_SMITH_TENSOR_H
#define __SRC_SMITH_TENSOR_H

#include <stddef.h>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <unordered_set>
#include <src/util/math/matrix.h>
#include <src/util/math/matop.h>
#include <src/util/prim_op.h>
#include <src/smith/storage.h>
#include <src/smith/indexrange.h>
#include <src/smith/loopgenerator.h>

namespace bagel {
namespace SMITH {

template <typename DataType>
class Tensor_ {
  public:
    using MatType = typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type;
    using VecType = typename std::conditional<std::is_same<DataType,double>::value, VectorB, ZVectorB>::type;
  protected:
    std::vector<IndexRange> range_;
    std::shared_ptr<Storage<DataType>> data_;
    int rank_;
    std::unordered_set<size_t> sparse_;

    mutable bool initialized_;

    bool allocated_;

  public:
    Tensor_() { }
    Tensor_(std::vector<IndexRange> in, const bool kramers = false, const std::unordered_set<size_t> sparse = std::unordered_set<size_t>(), const bool alloc = false);
    virtual ~Tensor_() { }

    Tensor_<DataType>& operator=(const Tensor_<DataType>& o);

    std::shared_ptr<Tensor_<DataType>> clone() const;
    std::shared_ptr<Tensor_<DataType>> copy() const;

    virtual void init() const { initialized_ = true; }

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << range_ << data_ << rank_ << sparse_ << initialized_ << allocated_;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      bool do_init;
      ar >> range_ >> data_ >> rank_ >> sparse_ >> do_init >> allocated_;
      initialized_ = false;
      if (do_init) init();
      if (allocated_ != data_->initialized())
        throw std::runtime_error("Allocation error when trying to load a serialized Smith tensor.");
    }

  public:
    void ax_plus_y(const DataType& a, const Tensor_<DataType>& o) { data_->ax_plus_y(a, o.data_); }
    void ax_plus_y(const DataType& a, std::shared_ptr<const Tensor_<DataType>> o) { ax_plus_y(a, *o); }

    void scale(const DataType& a) { data_->scale(a); }

    bool allocated() const { return allocated_; }
    void allocate();
    void fence() const { data_->fence(); }

    template<typename ...args>
    bool is_local(args&& ...p) const { return data_->is_local(std::forward<args>(p)...); }
    template<typename ...args>
    bool exists(args&& ...p) const { return sparse_.empty() || sparse_.count(generate_hash_key(std::forward<args>(p)...)); }

    DataType dot_product(const Tensor_<DataType>& o) const { return data_->dot_product(*o.data_); }
    DataType dot_product(std::shared_ptr<const Tensor_<DataType>> o) const { return dot_product(*o); }

    int rank() const { return rank_; }
    size_t size_alloc() const;

    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    double rms() const { return std::sqrt(detail::real(dot_product(*this))/size_alloc()); }

    std::vector<IndexRange> indexrange() const { return range_; }

    template<typename ...args>
    std::unique_ptr<DataType[]> get_block(args&& ...p) const {
      return data_->get_block(std::forward<args>(p)...);
    }

    template<typename ...args>
    void put_block(std::unique_ptr<DataType[]>& o, args&& ...p) {
      data_->put_block(o, std::forward<args>(p)...);
    }

    template<typename ...args>
    void add_block(std::unique_ptr<DataType[]>& o, args&& ...p) {
      data_->add_block(o, std::forward<args>(p)...);
    }

    template<typename ...args>
    size_t get_size(args&& ...p) const {
      return data_->blocksize(std::forward<args>(p)...);
    }

    void zero() { data_->zero(); }

    double orthog(const std::list<std::shared_ptr<const Tensor_<DataType>>> o);

    std::vector<DataType> diag() const;
    std::shared_ptr<MatType> matrix() const;
    std::shared_ptr<VecType> vectorb() const;

    // for Kramers tensors (does not do anything for standard tensors)
    void set_perm(const std::map<std::vector<int>, std::pair<double,bool>>& p) { data_->set_perm(p); }
    void set_stored_sectors(const std::list<std::vector<bool>>& s) { data_->set_stored_sectors(s); }
};

extern template class Tensor_<double>;
extern template class Tensor_<std::complex<double>>;

namespace CASPT2 { using Tensor = Tensor_<double>; }
namespace CASA { using Tensor = Tensor_<double>; }
namespace SPCASPT2 { using Tensor = Tensor_<double>; }
namespace MSCASPT2 { using Tensor = Tensor_<double>; }
namespace MRCI   { using Tensor = Tensor_<double>; }
namespace RelCASPT2 { using Tensor = Tensor_<std::complex<double>>; }
namespace RelCASA { using Tensor = Tensor_<std::complex<double>>; }
namespace RelSPCASPT2 { using Tensor = Tensor_<std::complex<double>>; }
namespace RelMSCASPT2 { using Tensor = Tensor_<std::complex<double>>; }
namespace RelMRCI   { using Tensor = Tensor_<std::complex<double>>; }

}
}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::SMITH::Tensor_<double>)
BOOST_CLASS_EXPORT_KEY(bagel::SMITH::Tensor_<std::complex<double>>)

#endif
