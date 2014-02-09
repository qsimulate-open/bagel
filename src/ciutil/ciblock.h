//
// BAGEL - Parallel electron correlation program.
// Filename: ciblock.h
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
//

#ifndef __SRC_CIUTIL_CIBLOCK_H
#define __SRC_CIUTIL_CIBLOCK_H

#include <src/ciutil/bitutil.h>
#include <src/ciutil/cistring.h>

namespace bagel {

template<class StringType>
class CIBlockInfo {
  protected:
    std::shared_ptr<const StringType> astrings_;
    std::shared_ptr<const StringType> bstrings_;

  public:
    CIBlockInfo() { }
    CIBlockInfo(std::shared_ptr<const StringType> ast, std::shared_ptr<const StringType> bst)
      : astrings_(ast), bstrings_(bst) {
      static_assert(std::is_base_of<CIString_base, StringType>::value, "illegal StringType specified");
      assert(astrings_->norb() == bstrings_->norb());
    }
    virtual ~CIBlockInfo() { }

    size_t size() const { return lena() * lenb(); }
    size_t lena() const { return astrings_->size(); }
    size_t lenb() const { return bstrings_->size(); }
    size_t norb() const { return astrings_->norb(); }
    int nspin() const { return nelea() - neleb(); }
    int nelea() const { return astrings_->nele(); }
    int neleb() const { return bstrings_->nele(); }

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    template <int spin>
    size_t lexical(const std::bitset<nbit__>& bit) const {
      return spin == Alpha ? astrings_->lexical(bit) : bstrings_->lexical(bit);
    }

    size_t index(const std::bitset<nbit__>& bbit, const std::bitset<nbit__>& abit) const {
      return bstrings_->lexical_zero(bbit) + astrings_->lexical_zero(abit) * lenb();
    }

    const std::shared_ptr<const StringType>& stringsa() const { return astrings_; }
    const std::shared_ptr<const StringType>& stringsb() const { return bstrings_; }

    template<int spin>
    const std::shared_ptr<const StringType>& strings() const { return spin == 0 ? astrings_ : bstrings_; }

    const std::bitset<nbit__>& string_bits_a(int i) const { return astrings_->strings(i); }
    const std::bitset<nbit__>& string_bits_b(int i) const { return bstrings_->strings(i); }
    const std::vector<std::bitset<nbit__>>& string_bits_a() const { return astrings_->strings(); }
    const std::vector<std::bitset<nbit__>>& string_bits_b() const { return bstrings_->strings(); }

    template<int spin>
    int sign(const std::bitset<nbit__>& bit, int i) const {
      return bagel::sign(bit, i)*(1-(((spin*nelea())&1)<<1));
    }
    static int sign(const std::bitset<nbit__>& bit, int i, int j) { return bagel::sign(bit, i, j); }

};


// Contains all the information for (a sub block of) the CI coefficient matrix
// but does NOT own the data
template <typename DataType, class StringType>
class CIBlock : public CIBlockInfo<StringType> {
  protected:
    DataType* data_ptr_;
    size_t offset_;

  public:
    CIBlock(std::shared_ptr<const StringType> astrings, std::shared_ptr<const StringType> bstrings, DataType* const data_ptr, const size_t o) :
      CIBlockInfo<StringType>(astrings, bstrings), data_ptr_(data_ptr), offset_(o) {
    }
    virtual ~CIBlock() { }

    DataType* data() { return data_ptr_; }
    const DataType* data() const { return data_ptr_; }

    DataType& element(const size_t i) { return data_ptr_[i]; }
    const DataType& element(const size_t i) const { return data_ptr_[i]; }

    DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) { return element( this->index(bstring, astring) ); }
    const DataType& element(const std::bitset<nbit__> bstring, const std::bitset<nbit__> astring) const { return element( this->index(bstring, astring) ); }

};


template <typename DataType, class StringType>
class CIBlock_alloc : public CIBlock<DataType, StringType> {
  protected:
    std::unique_ptr<DataType[]> data_;
  public:
    CIBlock_alloc(std::shared_ptr<const StringType> astrings, std::shared_ptr<const StringType> bstrings) : CIBlock<DataType,StringType>(astrings, bstrings, 0, 0) {
      data_ = std::unique_ptr<DataType[]>(new DataType[this->size()]);
      std::fill_n(data_.get(), this->size(), DataType(0.0));
      this->data_ptr_ = data_.get();
    }
};

}

#endif
