//
// BAGEL - Parallel electron correlation program.
// Filename: data.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_PERIODIC_PDATA_H
#define __SRC_PERIODIC_PDATA_H

#include <src/math/zmatrix.h>

namespace bagel {

/* Store data in direct space */
class PData {
  protected:
    int blocksize_;
    int nblock_;

    std::vector<std::shared_ptr<ZMatrix>> pdata_;     // (g, i, j)

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & blocksize_ & nblock_ & pdata_;
    }

  public:
    PData() { }
    PData(const int bsize, const int nblock) : blocksize_(bsize), nblock_(nblock) {
      pdata_.resize(nblock);
      for (int i = 0; i != nblock; ++i) {
        auto block = std::make_shared<ZMatrix>(bsize, bsize);
        block->zero();
        pdata_[i] = block;
      }
    }

    ~PData() { }

    const int blocksize() const { return blocksize_; }
    const int nblock() const { return nblock_; }

    std::shared_ptr<ZMatrix> operator[] (int i) { assert(i < nblock_ && i >= 0); return pdata_[i]; };

    std::vector<std::shared_ptr<ZMatrix>> pdata() const { return pdata_; }
    std::shared_ptr<ZMatrix> pdata(const int i) const { return pdata_[i]; }

    void zero() {
      for (auto& block : pdata_) block->zero();
    }

    void allreduce() {
      for (auto& block : pdata_) block->allreduce();
    }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::PData)

#endif
