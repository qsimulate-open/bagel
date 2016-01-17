//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/civector.cc
// Copyright (C) 2013 Toru Shiozaki
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


#include <iomanip>
#include <unordered_map>

#include <src/ci/ras/civector.h>

namespace bagel {
namespace RAS {
  /** \brief Implementation for application of \f$\hat S^2\f$

      \detail Implemented as \f[\hat S^2 = \hat S_z^2 + \hat S_z + \hat S_-\hat S_+\f]
              where \f[ \hat S_-\hat S_+ = n_\beta - \sum_{i,j} j_\alpha^\dagger i_\alpha i_\beta^\dagger j_\beta \f] */
  void spin_impl(const RASCivecView cc, RASCivecView out);

  /// Spin lowering operator \f$ \hat S_- = \sum_i i^\dagger_\beta i_\alpha \f$
  void spin_lower_impl(const RASCivecView cc, RASCivecView out);

  /// Spin raising operator \f$ \hat S_+ = \sum_i i^\dagger_\alpha i_\beta \f$
  void spin_raise_impl(const RASCivecView cc, RASCivecView out);
}
}
