//
// BAGEL - Parallel electron correlation program.
// Filename: phi_ijk_lists.h
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

#ifndef __SRC_ASD_DMRG_PHI_IJK_LISTS_H
#define __SRC_ASD_DMRG_PHI_IJK_LISTS_H

#include <src/ciutil/citraits.h>
#include <src/ras/civector.h>

namespace bagel {

/// A list of excitation lists of either \f$|I'\rangle = (-1)^\phi \hat{i}^\dagger\hat{j}^\dagger\hat{k}|I\rangle\f$
/// or its conjugate.
class PhiIJKLists {
  public: struct PhiIJK {
    size_t source;
    int sign;
    int i;
    int j;
    int k;
    PhiIJK(size_t s, int si, int _i, int _j, int _k) : source(s), sign(si), i(_i), j(_j), k(_k) {}
  };

  protected:
    std::vector<std::vector<PhiIJK>> data_;

  public:
/** Constructor that builds either (i) \f$|I'\rangle = (-1)^\phi \hat{i}^\dagger\hat{j}^\dagger\hat{k}|I\rangle\f$
 *  or (ii) \f$|I'\rangle = (-1)^\phi \hat{i}^\dagger\hat{j}\hat{k}|I\rangle\f$
    @param conjugate true does (i), false does (ii)
 */
    PhiIJKLists(std::shared_ptr<const CIStringSet<RASString>> source_stringspace,
                std::shared_ptr<const CIStringSet<RASString>> target_stringspace, const bool conjugate);

    const std::vector<PhiIJK>& data(const int ia) const { return data_.at(ia); }
};

}

#endif
