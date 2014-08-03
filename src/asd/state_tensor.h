//
// BAGEL - Parallel electron correlation program.
// Filename: state_tensor.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
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


#ifndef __SRC_ASD_STATE_TENSOR_H
#define __SRC_ASD_STATE_TENSOR_H

#include <src/math/matrix.h>
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

};

}

#endif
