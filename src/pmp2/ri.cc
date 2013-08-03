//
// BAGEL - Parallel electron correlation program.
// Filename: ri.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <src/pmp2/pmp2.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/ptildex.h>
#include <src/pscf/phcore.h>
#include <src/pscf/pcoeff.h>
#include <src/pscf/pcompcabsfile.h>

using namespace std;
using namespace bagel;

#ifdef HAVE_LIBSLATER

// create RI basis

typedef shared_ptr<PMatrix1e> RefMatrix;
typedef shared_ptr<PGeometry> RefGeom;
typedef shared_ptr<PHcore> RefHcore;
typedef shared_ptr<PCoeff> RefCoeff;

std::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> PMP2::generate_RI() {

  // Form RI space which is a union of OBS and CABS.
  RefGeom newgeom(new PGeometry(*geom_));
  union_geom_ = newgeom;
  union_geom_->merge_obs_aux();

  shared_ptr<POverlap> union_overlap(new POverlap(union_geom_));
  shared_ptr<PTildeX> ri_coeff(new PTildeX(union_overlap));

  ri_entire_ = ri_coeff;

  RefMatrix ri_o(new PMatrix1e(ri_entire_, make_pair(0, geom_->nbasis())));
  RefMatrix ri_c(new PMatrix1e(ri_entire_, make_pair(geom_->nbasis(), geom_->nbasis()+geom_->naux())));

  pair<RefMatrix, RefMatrix> ri_o_pair = ri_o->split(geom_->nbasis(), geom_->naux());
  pair<RefMatrix, RefMatrix> ri_c_pair = ri_c->split(geom_->nbasis(), geom_->naux());

  return make_tuple(ri_o_pair.first, ri_o_pair.second, ri_c_pair.first, ri_c_pair.second);
}


std::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> PMP2::generate_fock_weighted_RI() const {

  RefMatrix barefock(new PMatrix1e(*ri_entire_ % (*ao_hJ_-*ao_K_) * *ri_entire_));

  RefMatrix fock(new PMatrix1e(*ri_entire_ * *barefock));

  RefMatrix fock_o(new PMatrix1e(fock, make_pair(0, geom_->nbasis())));
  RefMatrix fock_c(new PMatrix1e(fock, make_pair(geom_->nbasis(), geom_->nbasis()+geom_->naux())));

  pair<RefMatrix, RefMatrix> fock_o_pair = fock_o->split(geom_->nbasis(), geom_->naux());
  pair<RefMatrix, RefMatrix> fock_c_pair = fock_c->split(geom_->nbasis(), geom_->naux());

  return make_tuple(fock_o_pair.first, fock_o_pair.second, fock_c_pair.first, fock_c_pair.second);
}

#endif
