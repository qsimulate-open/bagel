//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/civector_impl.h
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


#ifndef __BAGEL_RAS_RASCIVECTOR_IMPL_H
#define __BAGEL_RAS_RASCIVECTOR_IMPL_H

#include <list>
#include <src/util/math/algo.h>
#include <src/ci/ras/civector_base.h>
#include <src/ci/fci/dvector_base.h>

namespace bagel {

// partial specialization of CIBlock (ciutil/ciblock.h)
template<typename DataType>
using RASBlock = CIBlock<DataType, RASString>;
template<typename DataType>
using RASBlock_alloc = CIBlock_alloc<DataType, RASString>;

template <typename DataType>
class RASCivector_impl : public RASCivector_base<RASBlock<DataType>> {
  public:
    using DetType = RASDeterminants;
    using RBlock = RASBlock<DataType>;
    using LocalizedType = std::true_type;
    using RASCivector_base<RASBlock<DataType>>::block;

  protected:
    using RASCivector_base<RASBlock<DataType>>::blocks_;
    using RASCivector_base<RASBlock<DataType>>::det_;

    virtual std::shared_ptr<RASCivector_impl<DataType>> spin_() const = 0;

    template<typename T>
    void transpose_impl(std::shared_ptr<T>& out, std::shared_ptr<const RASDeterminants> det = nullptr) const {
      if (!det) det = det_->transpose();
      const int phase = 1 - (((det->nelea()*det->neleb())%2) << 1);
      out = std::make_shared<T>(det);
      this->for_each_block([&out, &phase] (std::shared_ptr<const RBlock> b) {
        blas::transpose(b->data(), b->lenb(), b->lena(), out->block(b->stringsa(), b->stringsb())->data(), static_cast<double>(phase));
      });
    }

  public:
    RASCivector_impl(std::shared_ptr<const RASDeterminants> det) : RASCivector_base<RBlock>(det) {}

    virtual DataType* data() = 0;
    virtual const DataType* data() const = 0;

    // Copy assignment
    RASCivector_impl<DataType>& operator=(const RASCivector_impl<DataType>& o) {
      assert(*o.det_ == *det_);
      std::copy_n(o.data(), size(), data());
      return *this;
    }

    // Element-wise access. Beware: very slow!
    DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) {
      return block(bstring, astring)->element(bstring, astring);
    }
    const DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) const {
      return block(bstring, astring)->element(bstring, astring);
    }

    size_t size() const { return det_->size(); }
    void fill(const double a) { std::fill_n(data(), size(), a); }
    void zero() { fill(0.0); }

    // Safe for any structure of blocks.
    template<typename T>
    DataType dot_product(const T& o) const {
      assert(det_->nelea() == o.det()->nelea() && det_->neleb() == o.det()->neleb() && det_->norb() == o.det()->norb());
      DataType out(0.0);
      this->for_each_block( [&out, &o] (std::shared_ptr<const RBlock> b) {
        std::shared_ptr<const RBlock> j = o.block(b->stringsb(), b->stringsa());
        if (j) out += blas::dot_product(b->data(), b->lena()*b->lenb(), j->data());
      } );
      return out;
    }
    template<typename T>
    DataType dot_product(std::shared_ptr<const T> o) const { return dot_product(*o); }

    DataType spin_expectation() const { return dot_product(*spin_()); }
    void spin_decontaminate(const double thresh = 1.0e-8);

    template<typename T>
    void ax_plus_y(const double a, const T& o) { blas::ax_plus_y_n(a, o.data(), size(), data()); }
    template<typename T>
    void ax_plus_y(const double a, std::shared_ptr<const T> o) { ax_plus_y(a, *o); }

    void scale(const DataType a) { blas::scale_n(a, data(), size()); }

    template<typename T>
    void project_out(const std::shared_ptr<const T>& o) { ax_plus_y(-dot_product(*o), *o); }

    double norm() const { return std::sqrt(blas::dot_product(data(), size(), data())); }
    double variance() const { return blas::dot_product(data(), size(), data())/size(); }
    double rms() const { return std::sqrt(variance()); }

    template<typename T>
    double orthog(std::list<std::shared_ptr<const T>> c) {
      for (auto& iter : c)
        project_out(iter);
      return normalize();
    }

    template<typename T>
    double orthog(const std::shared_ptr<const T>& o) {
      return orthog(std::list<std::shared_ptr<const T>>{o});
    }

    double normalize();
    void print(const double thr = 0.05) const;
    void synchronize(const int root = 0);
};

extern template class RASCivector_impl<double>;

}

#endif
