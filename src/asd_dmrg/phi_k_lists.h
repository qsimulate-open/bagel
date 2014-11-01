//
// BAGEL - Parallel electron correlation program.
// Filename: phi_k_lists.h
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

#ifndef __SRC_ASD_DMRG_PHI_K_LISTS_H
#define __SRC_ASD_DMRG_PHI_K_LISTS_H

#include <src/ciutil/citraits.h>
#include <src/ras/civector.h>

namespace bagel {

/// A set of initially blank SparseMatrix objects to help build matrices of the form \f$F(\beta',\beta) = f(i,j,\phi)\f$ where
/// \f$|\beta'\rangle = (-1)^\phi\hat E_{ij} |\beta\rangle\f$. Not restricted to \f$\beta\f$ strings.
class PhiKLists {
  public: struct PhiK {
    size_t source;
    size_t target;
    int sign;
    const RASString* target_space;
    PhiK(size_t s, size_t t, int sn, const RASString* ts) : source(s), target(t), sign(sn), target_space(ts) {}
  };

  protected:
    /// map<[orbital number], map<[space tag], vector<PhiK>>
    std::unordered_map<int, std::unordered_map<int, std::vector<PhiK>>> data_;

  public:
    PhiKLists(std::shared_ptr<const CIStringSet<RASString>> source_stringspace, std::shared_ptr<const CIStringSet<RASString>> target_stringspace);

    const std::unordered_map<int, std::vector<PhiK>>& data(const int k) const { return data_.at(k); }
};

}

#endif
