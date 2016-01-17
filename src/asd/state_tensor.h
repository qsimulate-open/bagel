//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: state_tensor.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
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


#ifndef __SRC_ASD_STATE_TENSOR_H
#define __SRC_ASD_STATE_TENSOR_H

#include <src/util/math/matrix.h>
#include <src/asd/dimersubspace.h>

namespace bagel {

class StateTensor {
  protected:
    using SparseMap = std::map<std::tuple<int, MonomerKey, MonomerKey>, btas::PTensor2<const double>>;

    // original data
    std::shared_ptr<const Matrix> data_;

    SparseMap sparse_;
    int norb_;
    int nstates_;

  public:
    StateTensor(std::shared_ptr<const Matrix> diab, const std::vector<DimerSubspace_base>& subspaces) : data_(diab), nstates_(diab->mdim()) {
      for (int istate = 0; istate != nstates_; ++istate) {
        for (auto& i : subspaces) {
          MonomerKey k0 = i.monomerkey<0>();
          MonomerKey k1 = i.monomerkey<1>();
          const int n0 = k0.nstates();
          const int n1 = k1.nstates();
          if (n0*n1 == 0) continue;

          const int offset = i.offset();
          PreAllocArray<const double> storage(diab->element_ptr(offset, istate), n0*n1);
          btas::CRange<2> range(n0, n1);
          btas::PTensor2<const double> block(range, storage);

          sparse_.emplace(std::make_tuple(istate, k0, k1), block);
        }
      }
    }

    void print() const {
      for (auto& i : sparse_) {
        std::string label = std::get<1>(i.first).to_string() + " " + std::get<2>(i.first).to_string();
        btas::print(i.second, label + " (state " + lexical_cast<std::string>(std::get<0>(i.first)) + ")");
        std::cout << std::endl;
      }
    }

    auto begin() -> decltype(sparse_.begin()) { return sparse_.begin(); }
    auto end() -> decltype(sparse_.end()) { return sparse_.end(); }
    auto begin() const -> decltype(sparse_.cbegin()) { return sparse_.cbegin(); }
    auto end() const -> decltype(sparse_.cend()) { return sparse_.cend(); }
    auto cbegin() const -> decltype(sparse_.cbegin()) { return sparse_.cbegin(); }
    auto cend() const -> decltype(sparse_.cend()) { return sparse_.cend(); }

    int nblocks() const { return sparse_.size(); }

    bool exist(const std::tuple<int, MonomerKey, MonomerKey>& i) const { return sparse_.find(i) != sparse_.end(); }

    btas::PTensor2<const double> get_block(const MonomerKey& i, const MonomerKey& j, const int k) const {
      return sparse_.at(std::make_tuple(k, i, j));
    }

    std::shared_ptr<Matrix> contract_statetensor(const std::array<MonomerKey,4>& keys, const int istate) const {
      auto& A  = keys[0]; auto& B  = keys[1];
      auto& Ap = keys[2]; auto& Bp = keys[3];

      assert(A == Ap);
      assert(exist(std::make_tuple(istate,A,B)));
      assert(exist(std::make_tuple(istate,Ap,Bp)));

      auto tensor = std::make_shared<btas::Tensor2<double>>(B.nstates(), Bp.nstates());
      btas::contract(1.0, sparse_.at(std::make_tuple(istate,A,B)), {0,1}, sparse_.at(std::make_tuple(istate,Ap,Bp)), {0,2}, 0.0, *tensor, {1,2});
      btas::CRange<2> range(tensor->extent(0)*tensor->extent(1), 1);
      MatView view(btas::make_view(range, tensor->storage()), /*localized*/true);

      return std::make_shared<Matrix>(view);
    }

};

}

#endif
