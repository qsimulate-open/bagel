//
// BAGEL - Parallel electron correlation program.
// Filename: ras/civector.cc
// Copyright (C) 2013 Toru Shiozaki
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


#include <iomanip>
#include <unordered_map>

#include <src/ras/civector.h>

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
