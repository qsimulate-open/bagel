//
// BAGEL - Parallel electron correlation program.
// Filename: alpha.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Matthew Kelley <matthewkelley2017@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#ifndef __SRC_REL_ALPHA_H
#define __SRC_REL_ALPHA_H

#include <src/util/zmatrix.h>
#include <src/rel/dfdata.h>

namespace bagel {

class Sigma {
  protected:
    std::shared_ptr<ZMatrix> data_;
  public:
    Sigma(const int i) : data_(new ZMatrix(2,2,true)) {
      if (i == DFData::Comp::L) {
        data_->element(0) = data_->element(3) = 1.0;
        data_->element(1) = data_->element(2) = 0.0;
      } else if (i == DFData::Comp::X) {
        data_->element(0) = data_->element(3) = 0.0;
        data_->element(1) = data_->element(2) = 1.0;
      } else if (i == DFData::Comp::Y) {
        data_->element(0) = data_->element(3) = 0.0;
        data_->element(1) = std::complex<int>(0,1);
        data_->element(2) = std::complex<int>(0,-1);
      } else if (i == DFData::Comp::Z) {
        data_->element(0) = 1.0;
        data_->element(3) = -1.0;
        data_->element(1) = data_->element(2) = 0.0;
      } else {
        assert(false);
      }
    }

    std::shared_ptr<const ZMatrix> data() const { return data_; }
};

class Alpha {
  protected:
    std::shared_ptr<ZMatrix> data_;
  public:
    Alpha(const int i) : data_(new ZMatrix(4,4,true)) {
      Sigma s(i);
      if (i == DFData::Comp::L) {
        data_->copy_block(0,0,2,2,s.data());
        data_->copy_block(2,2,2,2,s.data());
      else {
        data_->copy_block(2,0,2,2,s.data());
        data_->copy_block(0,2,2,2,s.data());
      }
      assert((*data_*data_).trace() == std::complex<double>(4.0));
    }

    std::shared_ptr<const ZMatrix> data() const { return data_; }

};

}

#endif
