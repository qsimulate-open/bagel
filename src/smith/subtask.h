//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: subtask.h
// Copyright (C) 2012 Toru Shiozaki
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
template<int N, int M, typename DataType>
class SubTask_ {
  protected:
    const std::array<const Index, N> block_index_;
    const std::array<std::shared_ptr<const Tensor_<DataType>>, M> in_;
    const std::shared_ptr<Tensor_<DataType>> out_;

  public:
    SubTask_(const std::array<const Index, N>& i, const std::array<std::shared_ptr<const Tensor_<DataType>>, M>& j, std::shared_ptr<Tensor_<DataType>>& k)
      : block_index_(i), in_(j), out_(k) { }

    virtual void compute() = 0;

    const Index& block(const size_t& i) const { return block_index_[i]; }
    const std::array<std::shared_ptr<const Tensor_<DataType>>, M>& in_tensors() const { return in_; }
    std::shared_ptr<const Tensor_<DataType>> in_tensor(const size_t& i) const { return in_[i]; }
    std::shared_ptr<Tensor_<DataType>> out_tensor() { return out_; }
};

// When we have multiple outputs
// N is the number of target indices, M is the number of input tensors, and L is the
// number of output tensors. The others are same.
template<int N, int M, int L, typename DataType>
class SubTask_Merged_ {
  protected:
    const std::array<const Index, N> block_index_;
    const std::array<std::shared_ptr<const Tensor_<DataType>>, M> in_;
    const std::array<std::shared_ptr<Tensor_<DataType>>, L> out_;

  public:
    SubTask_Merged_(const std::array<const Index, N>& i, const std::array<std::shared_ptr<const Tensor_<DataType>>, M>& j, const std::array<std::shared_ptr<Tensor_<DataType>>, L>& k)
      : block_index_(i), in_(j), out_(k) { }

    virtual void compute() = 0;

    const Index& block(const size_t& i) const { return block_index_[i]; }
    const std::array<std::shared_ptr<const Tensor_<DataType>>, M>& in_tensors() const { return in_; }
    const std::array<std::shared_ptr<Tensor_<DataType>>, L>& out_tensors() const { return out_; }
    std::shared_ptr<const Tensor_<DataType>> in_tensor(const size_t& i) const { return in_[i]; }
    std::shared_ptr<Tensor_<DataType>> out_tensor(const size_t& i) { return out_[i]; }
};

namespace CASPT2 { template<int N, int M> using SubTask = SubTask_<N,M,double>; }
namespace CASPT2 { template<int N, int M, int L> using SubTask_Merged = SubTask_Merged_<N,M,L,double>; }
namespace CASA { template<int N, int M> using SubTask = SubTask_<N,M,double>; }
namespace MSCASPT2 { template<int N, int M> using SubTask = SubTask_<N,M,double>; }
namespace MSCASPT2 { template<int N, int M, int L> using SubTask_Merged = SubTask_Merged_<N,M,L,double>; }
namespace SPCASPT2 { template<int N, int M> using SubTask = SubTask_<N,M,double>; }
namespace MRCI   { template<int N, int M> using SubTask = SubTask_<N,M,double>; }
namespace RelCASPT2 { template<int N, int M> using SubTask = SubTask_<N,M,std::complex<double>>; }
namespace RelCASA { template<int N, int M> using SubTask = SubTask_<N,M,std::complex<double>>; }
namespace RelMRCI   { template<int N, int M> using SubTask = SubTask_<N,M,std::complex<double>>; }

extern template class SubTask_Merged_<2,1,5,double>;
extern template class SubTask_Merged_<4,1,5,double>;
extern template class SubTask_Merged_<6,1,5,double>;
extern template class SubTask_Merged_<8,1,5,double>;
extern template class SubTask_Merged_<2,2,5,double>;
extern template class SubTask_Merged_<4,2,5,double>;
extern template class SubTask_Merged_<6,2,5,double>;
extern template class SubTask_Merged_<8,2,5,double>;
extern template class SubTask_<1,1,double>;
extern template class SubTask_<1,2,double>;
extern template class SubTask_<1,3,double>;
extern template class SubTask_<1,4,double>;
extern template class SubTask_<1,5,double>;
extern template class SubTask_<1,6,double>;
extern template class SubTask_<1,7,double>;
extern template class SubTask_<1,8,double>;
extern template class SubTask_<1,9,double>;
extern template class SubTask_<2,1,double>;
extern template class SubTask_<2,2,double>;
extern template class SubTask_<2,3,double>;
extern template class SubTask_<2,4,double>;
extern template class SubTask_<2,5,double>;
extern template class SubTask_<2,6,double>;
extern template class SubTask_<2,7,double>;
extern template class SubTask_<2,8,double>;
extern template class SubTask_<2,9,double>;
extern template class SubTask_<3,1,double>;
extern template class SubTask_<3,2,double>;
extern template class SubTask_<3,3,double>;
extern template class SubTask_<3,4,double>;
extern template class SubTask_<3,5,double>;
extern template class SubTask_<3,6,double>;
extern template class SubTask_<3,7,double>;
extern template class SubTask_<3,8,double>;
extern template class SubTask_<3,9,double>;
extern template class SubTask_<4,1,double>;
extern template class SubTask_<4,2,double>;
extern template class SubTask_<4,3,double>;
extern template class SubTask_<4,4,double>;
extern template class SubTask_<4,5,double>;
extern template class SubTask_<4,6,double>;
extern template class SubTask_<4,7,double>;
extern template class SubTask_<4,8,double>;
extern template class SubTask_<4,9,double>;
extern template class SubTask_<5,1,double>;
extern template class SubTask_<5,2,double>;
extern template class SubTask_<5,3,double>;
extern template class SubTask_<5,4,double>;
extern template class SubTask_<5,5,double>;
extern template class SubTask_<5,6,double>;
extern template class SubTask_<5,7,double>;
extern template class SubTask_<5,8,double>;
extern template class SubTask_<5,9,double>;
extern template class SubTask_<6,1,double>;
extern template class SubTask_<6,2,double>;
extern template class SubTask_<6,3,double>;
extern template class SubTask_<6,4,double>;
extern template class SubTask_<6,5,double>;
extern template class SubTask_<6,6,double>;
extern template class SubTask_<6,7,double>;
extern template class SubTask_<6,8,double>;
extern template class SubTask_<6,9,double>;
extern template class SubTask_<7,1,double>;
extern template class SubTask_<7,2,double>;
extern template class SubTask_<7,3,double>;
extern template class SubTask_<7,4,double>;
extern template class SubTask_<7,5,double>;
extern template class SubTask_<7,6,double>;
extern template class SubTask_<7,7,double>;
extern template class SubTask_<7,8,double>;
extern template class SubTask_<7,9,double>;
extern template class SubTask_<8,1,double>;
extern template class SubTask_<8,2,double>;
extern template class SubTask_<8,3,double>;
extern template class SubTask_<8,4,double>;
extern template class SubTask_<8,5,double>;
extern template class SubTask_<8,6,double>;
extern template class SubTask_<8,7,double>;
extern template class SubTask_<8,8,double>;
extern template class SubTask_<8,9,double>;
extern template class SubTask_<9,1,double>;
extern template class SubTask_<9,2,double>;
extern template class SubTask_<9,3,double>;
extern template class SubTask_<9,4,double>;
extern template class SubTask_<9,5,double>;
extern template class SubTask_<9,6,double>;
extern template class SubTask_<9,7,double>;
extern template class SubTask_<9,8,double>;
extern template class SubTask_<9,9,double>;

extern template class SubTask_<1,1,std::complex<double>>;
extern template class SubTask_<1,2,std::complex<double>>;
extern template class SubTask_<1,3,std::complex<double>>;
extern template class SubTask_<1,4,std::complex<double>>;
extern template class SubTask_<1,5,std::complex<double>>;
extern template class SubTask_<1,6,std::complex<double>>;
extern template class SubTask_<1,7,std::complex<double>>;
extern template class SubTask_<1,8,std::complex<double>>;
extern template class SubTask_<1,9,std::complex<double>>;
extern template class SubTask_<2,1,std::complex<double>>;
extern template class SubTask_<2,2,std::complex<double>>;
extern template class SubTask_<2,3,std::complex<double>>;
extern template class SubTask_<2,4,std::complex<double>>;
extern template class SubTask_<2,5,std::complex<double>>;
extern template class SubTask_<2,6,std::complex<double>>;
extern template class SubTask_<2,7,std::complex<double>>;
extern template class SubTask_<2,8,std::complex<double>>;
extern template class SubTask_<2,9,std::complex<double>>;
extern template class SubTask_<3,1,std::complex<double>>;
extern template class SubTask_<3,2,std::complex<double>>;
extern template class SubTask_<3,3,std::complex<double>>;
extern template class SubTask_<3,4,std::complex<double>>;
extern template class SubTask_<3,5,std::complex<double>>;
extern template class SubTask_<3,6,std::complex<double>>;
extern template class SubTask_<3,7,std::complex<double>>;
extern template class SubTask_<3,8,std::complex<double>>;
extern template class SubTask_<3,9,std::complex<double>>;
extern template class SubTask_<4,1,std::complex<double>>;
extern template class SubTask_<4,2,std::complex<double>>;
extern template class SubTask_<4,3,std::complex<double>>;
extern template class SubTask_<4,4,std::complex<double>>;
extern template class SubTask_<4,5,std::complex<double>>;
extern template class SubTask_<4,6,std::complex<double>>;
extern template class SubTask_<4,7,std::complex<double>>;
extern template class SubTask_<4,8,std::complex<double>>;
extern template class SubTask_<4,9,std::complex<double>>;
extern template class SubTask_<5,1,std::complex<double>>;
extern template class SubTask_<5,2,std::complex<double>>;
extern template class SubTask_<5,3,std::complex<double>>;
extern template class SubTask_<5,4,std::complex<double>>;
extern template class SubTask_<5,5,std::complex<double>>;
extern template class SubTask_<5,6,std::complex<double>>;
extern template class SubTask_<5,7,std::complex<double>>;
extern template class SubTask_<5,8,std::complex<double>>;
extern template class SubTask_<5,9,std::complex<double>>;
extern template class SubTask_<6,1,std::complex<double>>;
extern template class SubTask_<6,2,std::complex<double>>;
extern template class SubTask_<6,3,std::complex<double>>;
extern template class SubTask_<6,4,std::complex<double>>;
extern template class SubTask_<6,5,std::complex<double>>;
extern template class SubTask_<6,6,std::complex<double>>;
extern template class SubTask_<6,7,std::complex<double>>;
extern template class SubTask_<6,8,std::complex<double>>;
extern template class SubTask_<6,9,std::complex<double>>;
extern template class SubTask_<7,1,std::complex<double>>;
extern template class SubTask_<7,2,std::complex<double>>;
extern template class SubTask_<7,3,std::complex<double>>;
extern template class SubTask_<7,4,std::complex<double>>;
extern template class SubTask_<7,5,std::complex<double>>;
extern template class SubTask_<7,6,std::complex<double>>;
extern template class SubTask_<7,7,std::complex<double>>;
extern template class SubTask_<7,8,std::complex<double>>;
extern template class SubTask_<7,9,std::complex<double>>;
extern template class SubTask_<8,1,std::complex<double>>;
extern template class SubTask_<8,2,std::complex<double>>;
extern template class SubTask_<8,3,std::complex<double>>;
extern template class SubTask_<8,4,std::complex<double>>;
extern template class SubTask_<8,5,std::complex<double>>;
extern template class SubTask_<8,6,std::complex<double>>;
extern template class SubTask_<8,7,std::complex<double>>;
extern template class SubTask_<8,8,std::complex<double>>;
extern template class SubTask_<8,9,std::complex<double>>;
extern template class SubTask_<9,1,std::complex<double>>;
extern template class SubTask_<9,2,std::complex<double>>;
extern template class SubTask_<9,3,std::complex<double>>;
extern template class SubTask_<9,4,std::complex<double>>;
extern template class SubTask_<9,5,std::complex<double>>;
extern template class SubTask_<9,6,std::complex<double>>;
extern template class SubTask_<9,7,std::complex<double>>;
extern template class SubTask_<9,8,std::complex<double>>;
extern template class SubTask_<9,9,std::complex<double>>;

}
}

#endif
