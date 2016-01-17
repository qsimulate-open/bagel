//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/sparse_ij.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __SRC_RAS_SPARSE_IJ_H
#define __SRC_RAS_SPARSE_IJ_H

#include <src/ci/ciutil/citraits.h>
#include <src/ci/ras/civector.h>
#include <src/util/math/sparsematrix.h>

namespace bagel {

/// A set of initially blank SparseMatrix objects to help build matrices of the form \f$F(\beta',\beta) = f(i,j,\phi)\f$ where
/// \f$|\beta'\rangle = (-1)^\phi\hat E_{ij} |\beta\rangle\f$. Not restricted to \f$\beta\f$ strings.
class Sparse_IJ {
  struct SparseIJKey {
    int i;
    int j;
    int sign;
    double* ptr;
    SparseIJKey(int ii, int jj, int s, double* p) : i(ii), j(jj), sign(s), ptr(p) {}
  };
  protected:
    std::map<std::pair<int, int>, std::tuple<std::shared_ptr<SparseMatrix>, std::vector<SparseIJKey>>> data_;

  public:
    Sparse_IJ(std::shared_ptr<const CIStringSet<RASString>> source_stringspace, std::shared_ptr<const CIStringSet<RASString>> target_stringspace);

    const std::tuple<std::shared_ptr<SparseMatrix>, std::vector<SparseIJKey>>& data(const int target_tag, const int source_tag) const { return data_.at({target_tag, source_tag}); }
    const std::shared_ptr<SparseMatrix>& sparse_matrix(const int target_tag, const int source_tag) const { return std::get<0>(data_.at({target_tag, source_tag})); }
    const std::vector<SparseIJKey>& sparse_data(const int target_tag, const int source_tag) const { return std::get<1>(data_.at({target_tag, source_tag})); }
};

}

#endif
