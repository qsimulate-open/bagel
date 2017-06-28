//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: storagekramers.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_SMITH_STORAGEKRAMERS_H
#define __SRC_SMITH_STORAGEKRAMERS_H

#include <list>
#include <sstream>
#include <src/util/kramers.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {

template<typename DataType>
class StorageKramers : public StorageIncore<DataType> {
  protected:
    std::list<std::vector<bool>> stored_sectors_;
    std::map<std::vector<int>, std::pair<double,bool>> perm_;

    template<int N>
    std::pair<std::vector<int>, std::pair<double,bool>> find_permutation(const KTag<N>& tag) const {
      auto trans = std::make_pair(std::vector<int>{0}, std::make_pair(0.0,false));
      for (auto& i : perm_)
        for (auto& j : stored_sectors_) {
          const KTag<N> ctag(j);
          if (tag == ctag.perm(i.first)) {
            trans = i;
            goto end;
          }
        }
      end:
      return trans;
    }

    template<typename... args>
    std::unique_ptr<DataType[]> get_block_(args&& ...key) const {
      constexpr const size_t N = sizeof...(key);
      std::vector<Index> indices = arg_convert(key...);
      std::vector<bool> kramers(N);
      for (int i = 0; i != N; ++i)
        kramers[i] = indices[i].kramers();

      // if this block is stored return immediately
      auto iter = std::find(stored_sectors_.begin(), stored_sectors_.end(), kramers);
      if (iter != stored_sectors_.end())
        return RMAWindow<DataType>::rma_get(generate_hash_key(key...));

      // if not, first find the right permutation
      const KTag<N> tag(kramers);
      std::pair<std::vector<int>, std::pair<double,bool>> trans = find_permutation(tag);
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
        const std::unique_ptr<DataType[]> data = RMAWindow<DataType>::rma_get(generate_hash_key(dindices));

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
    void put_block_(const std::unique_ptr<DataType[]>& dat, const args& ...key) {
      put_block_(dat, arg_convert(key...));
    }

    void put_block_(const std::unique_ptr<DataType[]>& dat, std::vector<Index> indices) {
#ifndef NDEBUG
      std::vector<bool> kramers(indices.size());
      for (int i = 0; i != indices.size(); ++i)
        kramers[i] = indices[i].kramers();
      if (std::find(stored_sectors_.begin(), stored_sectors_.end(), kramers) == stored_sectors_.end())
        throw std::logic_error("Kramers::put_block should only be called for existing blocks");
#endif
      RMAWindow<DataType>::rma_put(dat, generate_hash_key(indices));
    }

    template<typename... args>
    void add_block_(const std::unique_ptr<DataType[]>& dat, const args& ...key) {
      add_block_(dat, arg_convert(key...));
    }

    void add_block_(const std::unique_ptr<DataType[]>& dat, std::vector<Index> indices) {
#ifndef NDEBUG
      std::vector<bool> kramers(indices.size());
      for (int i = 0; i != indices.size(); ++i)
        kramers[i] = indices[i].kramers();
      if (std::find(stored_sectors_.begin(), stored_sectors_.end(), kramers) == stored_sectors_.end())
        throw std::logic_error("Kramers::add_block should only be called for existing blocks");
#endif
      RMAWindow<DataType>::rma_add(dat, generate_hash_key(indices));
    }

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<StorageIncore<DataType>>(*this) & stored_sectors_ & perm_;
    }

  public:
    StorageKramers() { }
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

    void put_block(const std::unique_ptr<DataType[]>& dat) override;
    void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0) override;
    void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) override;
    void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) override;
    void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) override;
    void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4) override;
    void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5) override;
    void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5, const Index& i6) override;
    void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5, const Index& i6, const Index& i7) override;
    void put_block(const std::unique_ptr<DataType[]>& dat, const std::vector<Index> indices);

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

    void set_perm(const std::map<std::vector<int>, std::pair<double,bool>>& p) override { perm_ = p; }
    void set_stored_sectors(const std::list<std::vector<bool>>& s) override { stored_sectors_ = s; }

};

extern template class StorageKramers<double>;
extern template class StorageKramers<std::complex<double>>;

template<typename DataType>
using StorageK = StorageKramers<DataType>;

}
}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::SMITH::StorageKramers<double>)
BOOST_CLASS_EXPORT_KEY(bagel::SMITH::StorageKramers<std::complex<double>>)

#endif
