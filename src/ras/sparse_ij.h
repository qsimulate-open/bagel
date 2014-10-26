//
// BAGEL - Parallel electron correlation program.
// Filename: ras/sparse_ij.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __SRC_RAS_SPARSE_IJ_H
#define __SRC_RAS_SPARSE_IJ_H

#include <src/ciutil/citraits.h>
#include <src/ras/civector.h>
#include <src/math/sparsematrix.h>

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
