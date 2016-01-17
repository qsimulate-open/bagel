//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: futuretensor.h
// Copyright (C) 2014 Toru Shiozaki
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
      TATensor<DataType,N>::operator=(std::move(*tensor_->template tiledarray<N>()));
    }
};

template<typename DataType, int N>
class FutureTATensor_new : public TATensor<DataType,N> {
  protected:
    std::shared_ptr<TATensor<DataType,N>> tensor_;
    mutable std::shared_ptr<Task> init_;

  public:
    FutureTATensor_new(std::shared_ptr<TATensor<DataType,N>> i, std::shared_ptr<Task> j) : TATensor<DataType,N>(*i), tensor_(i), init_(j) { }

    void init() override {
      this->initialized_ = true;
      init_->compute();
      TATensor<DataType,N>::operator=(std::move(*tensor_));
    }
};

}
}

#endif
