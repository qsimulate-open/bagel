//
// BAGEL - Parallel electron correlation program.
// Filename: subtask.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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

#ifndef __SRC_SMITH_SUBTASK_H
#define __SRC_SMITH_SUBTASK_H

namespace bagel {
namespace SMITH {

#include <array>
#include <memory>
#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>

// this is the class for single contraction

// N is the number of target indices, M is the number of input tensors.
// I am using template, because vectors are slow due to dynamic allocation, and
// subtask generation can be time consuming. By using array, the allocating will be
// at compile time
template<int N, int M, typename T>
class SubTask {
  protected:
    const std::array<const Index, N> block_index_;
    const std::array<std::shared_ptr<const Tensor<T>>, M> in_;
    const std::shared_ptr<Tensor<T>> out_;

  public:
    SubTask(const std::array<const Index, N>& i, const std::array<std::shared_ptr<const Tensor<T>>, M>& j, std::shared_ptr<Tensor<T>>& k)
      : block_index_(i), in_(j), out_(k) { }

    virtual void compute() = 0;

    const Index& block(const size_t& i) const { return block_index_[i]; }
    const std::shared_ptr<const Tensor<T>>& in_tensor(const size_t& i) const { return in_[i]; }
    const std::shared_ptr<Tensor<T>>& out_tensor() const { return out_; }
};

}
}

#endif
