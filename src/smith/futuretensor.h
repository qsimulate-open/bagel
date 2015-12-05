//
// BAGEL - Parallel electron correlation program.
// Filename: futuretensor.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __SRC_SMITH_FUTURETENSOR_H
#define __SRC_SMITH_FUTURETENSOR_H

#include <src/smith/tensor.h>
#include <src/smith/task.h>

namespace bagel {
namespace SMITH {

template<typename DataType, int N>
class FutureTATensor : public TATensor<DataType,N> {
  protected:
    using TATensor<DataType,N>::initialized_;

  protected:
    std::shared_ptr<Tensor_<DataType>> tensor_;
    mutable std::shared_ptr<Task> init_;

  public:
    FutureTATensor(const TATensor<DataType,N>& i,  std::shared_ptr<Tensor_<DataType>> t, std::shared_ptr<Task> j) : TATensor<DataType,N>(i), tensor_(t), init_(j) { }

    void init() override {
      init_->compute();
      // shallow copy
      TATensor<DataType,N>::operator=(*tensor_->template tiledarray<N>());
    }
};

}
}

#endif
