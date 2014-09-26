//
// BAGEL - Parallel electron correlation program.
// Filename: kdata.h
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


#ifndef __SRC_PERIODIC_KDATA_H
#define __SRC_PERIODIC_KDATA_H

#include <src/math/zmatrix.h>

namespace bagel {

/* Store data in reciprocal space */
class KData {
  protected:
    int blocksize_;
    int nblock_;

    std::vector<std::shared_ptr<ZMatrix>> kdata_;     // (k, i, j)

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & blocksize_ & nblock_ & kdata_;
    }

  public:
    KData() { }
    KData(const int bsize, const int nblock) : blocksize_(bsize), nblock_(nblock) {
      kdata_.resize(nblock);
      for (int i = 0; i != nblock; ++i) {
        auto block = std::make_shared<ZMatrix>(bsize, bsize);
        block->zero();
        kdata_[i] = block;
      }
    }

    ~KData() { }

    const int blocksize() const { return blocksize_; }
    const int nblock() const { return nblock_; }

    std::vector<std::shared_ptr<ZMatrix>> kdata() const { return kdata_; }
    std::shared_ptr<ZMatrix> kdata(const int i) const { return kdata_[i]; }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::KData)

#endif
