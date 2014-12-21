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

#ifndef __SRC_SMITH_SUBTASK_H
#define __SRC_SMITH_SUBTASK_H

#include <array>
#include <memory>
#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>

namespace bagel {
namespace SMITH {

// this is the class for single contraction

// N is the number of target indices, M is the number of input tensors.
// I am using template, because vectors are slow due to dynamic allocation, and
// subtask generation can be time consuming. By using array, the allocating will be
// at compile time
template<int N, int M>
class SubTask {
  protected:
    const std::array<const Index, N> block_index_;
    const std::array<std::shared_ptr<const Tensor>, M> in_;
    const std::shared_ptr<Tensor> out_;

  public:
    SubTask(const std::array<const Index, N>& i, const std::array<std::shared_ptr<const Tensor>, M>& j, std::shared_ptr<Tensor>& k)
      : block_index_(i), in_(j), out_(k) { }

    virtual void compute() = 0;

    const Index& block(const size_t& i) const { return block_index_[i]; }
    const std::shared_ptr<const Tensor>& in_tensor(const size_t& i) const { return in_[i]; }
    const std::shared_ptr<Tensor>& out_tensor() const { return out_; }
};

extern template class SubTask<1,1>;
extern template class SubTask<1,2>;
extern template class SubTask<1,3>;
extern template class SubTask<1,4>;
extern template class SubTask<1,5>;
extern template class SubTask<1,6>;
extern template class SubTask<1,7>;
extern template class SubTask<1,8>;
extern template class SubTask<1,9>;
extern template class SubTask<2,1>;
extern template class SubTask<2,2>;
extern template class SubTask<2,3>;
extern template class SubTask<2,4>;
extern template class SubTask<2,5>;
extern template class SubTask<2,6>;
extern template class SubTask<2,7>;
extern template class SubTask<2,8>;
extern template class SubTask<2,9>;
extern template class SubTask<3,1>;
extern template class SubTask<3,2>;
extern template class SubTask<3,3>;
extern template class SubTask<3,4>;
extern template class SubTask<3,5>;
extern template class SubTask<3,6>;
extern template class SubTask<3,7>;
extern template class SubTask<3,8>;
extern template class SubTask<3,9>;
extern template class SubTask<4,1>;
extern template class SubTask<4,2>;
extern template class SubTask<4,3>;
extern template class SubTask<4,4>;
extern template class SubTask<4,5>;
extern template class SubTask<4,6>;
extern template class SubTask<4,7>;
extern template class SubTask<4,8>;
extern template class SubTask<4,9>;
extern template class SubTask<5,1>;
extern template class SubTask<5,2>;
extern template class SubTask<5,3>;
extern template class SubTask<5,4>;
extern template class SubTask<5,5>;
extern template class SubTask<5,6>;
extern template class SubTask<5,7>;
extern template class SubTask<5,8>;
extern template class SubTask<5,9>;
extern template class SubTask<6,1>;
extern template class SubTask<6,2>;
extern template class SubTask<6,3>;
extern template class SubTask<6,4>;
extern template class SubTask<6,5>;
extern template class SubTask<6,6>;
extern template class SubTask<6,7>;
extern template class SubTask<6,8>;
extern template class SubTask<6,9>;
extern template class SubTask<7,1>;
extern template class SubTask<7,2>;
extern template class SubTask<7,3>;
extern template class SubTask<7,4>;
extern template class SubTask<7,5>;
extern template class SubTask<7,6>;
extern template class SubTask<7,7>;
extern template class SubTask<7,8>;
extern template class SubTask<7,9>;
extern template class SubTask<8,1>;
extern template class SubTask<8,2>;
extern template class SubTask<8,3>;
extern template class SubTask<8,4>;
extern template class SubTask<8,5>;
extern template class SubTask<8,6>;
extern template class SubTask<8,7>;
extern template class SubTask<8,8>;
extern template class SubTask<8,9>;
extern template class SubTask<9,1>;
extern template class SubTask<9,2>;
extern template class SubTask<9,3>;
extern template class SubTask<9,4>;
extern template class SubTask<9,5>;
extern template class SubTask<9,6>;
extern template class SubTask<9,7>;
extern template class SubTask<9,8>;
extern template class SubTask<9,9>;

}
}

#endif
