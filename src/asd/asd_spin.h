//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_spin.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifndef __asd_asd_spin_h
#define __asd_asd_spin_h

#include <src/util/math/sparsematrix.h>
#include <src/asd/coupling.h>
#include <src/asd/dimersubspace.h>

namespace bagel {

template <typename VecType>
class ASDSpinMap {
  protected:
    std::map<std::pair<int,int>, double> coords_;
  public:
    ASDSpinMap() { }

    const std::map<std::pair<int,int>, double>& coords() const { return coords_; }
    template<typename... Args>
    void emplace(Args&&... t) { coords_.emplace(std::forward<Args>(t)...); }

    void couple_blocks(const DimerSubspace<VecType>& AB, const DimerSubspace<VecType>& ApBp);
    void diagonal_block(const DimerSubspace<VecType>& AB);
};

#define ASD_HEADERS
#include <src/asd/asd_spin_coupling.hpp>
#undef  ASD_HEADERS


class ASDSpin : public SparseMatrix {
  protected:
    const int max_spin_;

  public:
    template<typename T>
    ASDSpin(const int dimension, const ASDSpinMap<T>& spinmap, const int max) : SparseMatrix(dimension, dimension, spinmap.coords()), max_spin_(max) {}

    int max() const { return max_spin_; }

    std::shared_ptr<Matrix> apply(const Matrix& o) const;
    void filter(Matrix& o, const int desired_spin) const;
};

}

#endif
