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

#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {

template<typename DataType>
class StorageKramers : public StorageIncore<DataType> {
  protected:
    template<typename... args>
    std::unique_ptr<DataType[]> get_block_(const args& ...key) const {
      // TODO implement
      return std::unique_ptr<DataType[]>();
    }

    template<typename... args>
    std::unique_ptr<DataType[]> move_block_(const args& ...key) {
      // TODO implement
      return std::unique_ptr<DataType[]>();
    }

    template<typename... args>
    void put_block_(std::unique_ptr<DataType[]>& dat, const args& ...key) {
      // TODO implement
    }

    template<typename... args>
    void add_block_(const std::unique_ptr<DataType[]>& dat, const args& ...key) {
      // TODO implement
    }

  public:
    // TODO temp constructor
    StorageKramers(const std::map<size_t, size_t>& size, bool init) : StorageIncore<DataType>(size, init) { }

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

    void ax_plus_y(const DataType& a, const StorageIncore<DataType>& o)                 { assert(false); }
    void ax_plus_y(const DataType& a, const std::shared_ptr<StorageIncore<DataType>> o) { assert(false); };
    DataType dot_product(const StorageIncore<DataType>& o) const                        { assert(false); return 0.0; }
};

extern template class StorageKramers<double>;
extern template class StorageKramers<std::complex<double>>;

template<typename DataType>
using StorageK = StorageKramers<DataType>;

}
}

#endif
