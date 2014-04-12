//
// BAGEL - Parallel electron correlation program.
// Filename: ecpbatch.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

#include <src/integral/ecp/ecpbatch.h>

using namespace bagel;
using namespace std;

void ECPBatch::integrate_angular() {//have to loop over atoms

  for (auto& aiter : mol_->atoms()) {
    shared_ptr<const ECP> aiter_ecp = aiter->ecp_parameters();
    RadialInt<AngularBatch, const shared_ptr<const ECP>, const array<shared_ptr<const Shell>,2>&>
                            radint(max_iter_, integral_thresh_, aiter_ecp, basisinfo_);
    *angulardata_ = radint.integral();
    ++angulardata_;
  }

}
