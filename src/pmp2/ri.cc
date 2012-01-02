/*
 * ri.cc
 *
 *  Created on: Nov 8, 2009
 *      Author: shiozaki
 */

#include <src/pmp2/pmp2.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/ptildex.h>
#include <src/pscf/phcore.h>
#include <src/pscf/pcoeff.h>
#include <src/util/pcompcabsfile.h>

using namespace std;

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

