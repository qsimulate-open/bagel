//
// BAGEL - Parallel electron correlation program.
// Filename: storagekramers.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_SMITH_STORAGEKRAMERS_H
#define __SRC_SMITH_STORAGEKRAMERS_H

#include <list>
#include <sstream>
#include <src/util/kramers.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {

namespace {
  void arg_convert_impl(std::vector<Index>& a) { }
  template<typename... args>
  void arg_convert_impl(std::vector<Index>& a, const Index& i, args... tail) {
    a.push_back(i);
    arg_convert_impl(a, tail...);
  }
  template<typename... args>
  std::vector<Index> arg_convert(args... p) {
    std::vector<Index> a;
    arg_convert_impl(a, p...);
    return a;
  }
}


template<typename DataType>
class StorageKramers : public StorageIncore<DataType> {
  protected:
    std::list<std::vector<bool>> stored_sector_;
    std::map<std::vector<int>, std::pair<double,bool>> perm_;

    template<typename... args>
    std::unique_ptr<DataType[]> get_block_(args&& ...key) const {
      constexpr const size_t N = sizeof...(key);
      std::vector<Index> indices = arg_convert(key...);
      std::vector<bool> kramers(N);
      for (int i = 0; i != N; ++i)
        kramers[i] = indices[i].kramers();

      // if this block is stored return immediately
      auto iter = std::find(stored_sector_.begin(), stored_sector_.end(), kramers);
      if (iter != stored_sector_.end())
        return StorageIncore<DataType>::get_block_(generate_hash_key(key...));

      // if not, first find the right permutation
      const KTag<N> tag(kramers);
      std::pair<std::vector<int>, std::pair<double,bool>>
        trans = std::make_pair(std::vector<int>{0}, std::make_pair(0.0,false));
      for (auto& i : perm_) {
        bool found = false;
        for (auto& j : stored_sector_) {
          const KTag<N> ctag(j);
          if (tag == ctag.perm(i.first)) {
            trans = i;
            found = true;
            break;
          }
        }
        if (found) break;
      }
      // allocate output area
      size_t buffersize = 1ull;
      for (auto& i : indices)
        buffersize *= i.size();
      std::unique_ptr<DataType[]> out(new DataType[buffersize]);

      if (trans.second.first != 0.0) {
        // reorder indices so that one can find the block. And retrieve.
        std::vector<Index> dindices(N);
        for (int i = 0; i != N; ++i)
          dindices[trans.first[i]] = indices[i];
        if (buffersize != this->blocksize(dindices)) {
          std::stringstream ss; ss << "incosistent : " << buffersize << " " << this->blocksize(dindices);
          throw std::logic_error(ss.str());
        }
        const std::unique_ptr<DataType[]> data = StorageIncore<DataType>::get_block_(generate_hash_key(dindices));

        // finally sort the date to the final format
        std::array<int,N> info, dim;
        for (int i = 0; i != N; ++i) {
          info[i] = trans.first[i];
          dim[i] = dindices[i].size();
        }
        sort_indices(info, trans.second.first, /*fac2*/0.0, data.get(), out.get(), dim);
        if (trans.second.second)
          blas::conj_n(out.get(), buffersize);
      } else {
        std::fill_n(out.get(), buffersize, 0.0);
      }
      return std::move(out);
    }

    template<typename... args>
    std::unique_ptr<DataType[]> move_block_(const args& ...key) {
      throw std::logic_error("StorageKramers is designed for input tensors");
      return std::unique_ptr<DataType[]>();
    }

    template<typename... args>
    void put_block_(std::unique_ptr<DataType[]>& dat, const args& ...key) {
      put_block_(dat, arg_convert(key...));
    }

    void put_block_(std::unique_ptr<DataType[]>& dat, std::vector<Index> indices) {
      StorageIncore<DataType>::put_block_(dat, generate_hash_key(indices));
      std::vector<bool> kramers(indices.size());
      for (int i = 0; i != indices.size(); ++i)
        kramers[i] = indices[i].kramers();
      if (std::find(stored_sector_.begin(), stored_sector_.end(), kramers) == stored_sector_.end())
        stored_sector_.push_back(kramers);
    }

    template<typename... args>
    void add_block_(const std::unique_ptr<DataType[]>& dat, const args& ...key) {
      throw std::logic_error("StorageKramers is designed for input tensors");
    }

  public:
    // TODO temp constructor
    StorageKramers(const std::map<size_t, size_t>& size, const bool init)
      : StorageIncore<DataType>(size, init) {
      for (auto& i : perm_)
        perm_.emplace(std::vector<int>(i.first.begin(), i.first.end()), i.second);
    }

    std::unique_ptr<DataType[]> get_block() const override;
    std::unique_ptr<DataType[]> get_block(const Index& i0) const override;
    std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1) const override;
    std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2) const override;
    std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3) const override;
    std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                          const Index& i4) const override;
    std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                          const Index& i4, const Index& i5) const override;
    std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                          const Index& i4, const Index& i5, const Index& i6) const override;
    std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                          const Index& i4, const Index& i5, const Index& i6, const Index& i7) const override;
    std::unique_ptr<DataType[]> get_block(std::vector<Index> i) const override { assert(false); return std::unique_ptr<DataType[]>(); }

    std::unique_ptr<DataType[]> move_block() override;
    std::unique_ptr<DataType[]> move_block(const Index& i0) override;
    std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1) override;
    std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2) override;
    std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3) override;
    std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                           const Index& i4) override;
    std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                           const Index& i4, const Index& i5) override;
    std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                           const Index& i4, const Index& i5, const Index& i6) override;
    std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                           const Index& i4, const Index& i5, const Index& i6, const Index& i7) override;

    void put_block(std::unique_ptr<DataType[]>& dat) override;
    void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0) override;
    void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) override;
    void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) override;
    void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) override;
    void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                     const Index& i4) override;
    void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                     const Index& i4, const Index& i5) override;
    void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                     const Index& i4, const Index& i5, const Index& i6) override;
    void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                     const Index& i4, const Index& i5, const Index& i6, const Index& i7) override;
    void put_block(std::unique_ptr<DataType[]>& dat, std::vector<Index> indices);

    void add_block(const std::unique_ptr<DataType[]>& dat) override;
    void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0) override;
    void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) override;
    void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) override;
    void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) override;
    void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4) override;
    void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5) override;
    void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5, const Index& i6) override;
    void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5, const Index& i6, const Index& i7) override;

    void ax_plus_y(const DataType& a, const StorageIncore<DataType>& o) {
      throw std::logic_error("StorageKramers::ax_plus_y illegal call");
    }
    DataType dot_product(const StorageIncore<DataType>& o) const {
      throw std::logic_error("StorageKramers::dot_product illegal call");
      return 0.0;
    }

    void set_perm(const std::map<std::vector<int>, std::pair<double,bool>>& p) override { perm_ = p; }
};

extern template class StorageKramers<double>;
extern template class StorageKramers<std::complex<double>>;

template<typename DataType>
using StorageK = StorageKramers<DataType>;

}
}

#endif
