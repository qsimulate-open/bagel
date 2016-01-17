//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ciblock.h
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
//

#ifndef __SRC_CIUTIL_CIBLOCK_H
#define __SRC_CIUTIL_CIBLOCK_H

#include <src/ci/ciutil/citraits.h>
#include <src/ci/ciutil/bitutil.h>

namespace bagel {

template<class StringType,
         class = typename std::enable_if<is_cistring<StringType>::value>::type
        >
class CIBlockInfo {
  public:
    typedef StringType string_type;

  protected:
    std::shared_ptr<const StringType> astrings_;
    std::shared_ptr<const StringType> bstrings_;

    size_t offset_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & astrings_ & bstrings_; }

  public:
    CIBlockInfo() { }
    CIBlockInfo(std::shared_ptr<const StringType> ast, std::shared_ptr<const StringType> bst, const size_t o = 0)
      : astrings_(ast), bstrings_(bst), offset_(o) {
      static_assert(std::is_base_of<CIString_base, StringType>::value, "illegal StringType specified");
      // norb should match or one of strings is dummy
      assert(astrings_->norb() == bstrings_->norb() || astrings_->norb()*bstrings_->norb() == 0);
    }
    virtual ~CIBlockInfo() { }

    virtual size_t size() const { return lena() * lenb(); }
    size_t lena() const { return astrings_->size(); }
    size_t lenb() const { return bstrings_->size(); }
    size_t norb() const { return astrings_->norb(); }
    int nspin() const { return nelea() - neleb(); }
    int nelea() const { return astrings_->nele(); }
    int neleb() const { return bstrings_->nele(); }
    size_t offset() const { return offset_; }

    bool empty() const { return size() == 0; }

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
  public:
    typedef DataType data_type;

  protected:
    DataType* data_ptr_;

  public:
    CIBlock(std::shared_ptr<const StringType> astrings, std::shared_ptr<const StringType> bstrings, DataType* const data_ptr, const size_t o) :
      CIBlockInfo<StringType>(astrings, bstrings, o), data_ptr_(data_ptr) {
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



// Contains and owns all the data and information for a sub block of the CI coefficient matrix
template <typename DataType, class StringType>
class DistCIBlock_alloc : public CIBlockInfo<StringType> {
  public:
    typedef DataType data_type;

  protected:
    const StaticDist dist_;

    // allocation size
    size_t astart_;
    size_t aend_;

    std::unique_ptr<DataType[]> local_;

    // Used during MPI routines
    const size_t block_offset_;

  public:
    DistCIBlock_alloc(std::shared_ptr<const StringType> astrings, std::shared_ptr<const StringType> bstrings, const size_t o) :
      CIBlockInfo<StringType>(astrings, bstrings), dist_(astrings->size(), mpi__->size()), block_offset_(o)
    {
      std::tie(astart_, aend_) = dist_.range(mpi__->rank());
      local_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      std::fill_n(local_.get(), size(), 0.0);
      mutex_ = std::vector<std::mutex>(asize());
    }
    // mutex for write accesses to local_
    mutable std::vector<std::mutex> mutex_;

    const StaticDist& dist() const { return dist_; }
    size_t block_offset() const { return block_offset_; }

    size_t asize() const { return aend_ - astart_; }
    size_t astart() const { return astart_; }
    size_t aend() const { return aend_; }

    size_t size() const override { return (aend_ - astart_) * this->lenb(); }
    size_t global_size() const { return this->lena() * this->lenb(); }

    DataType* local() { return local_.get(); }
    const DataType* local() const { return local_.get(); }
};

}

#endif
