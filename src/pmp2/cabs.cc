//
// BAGEL - Parallel electron correlation program.
// Filename: cabs.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/pmp2/pmp2.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/ptildex.h>
#include <src/pscf/phcore.h>
#include <src/pscf/pcoeff.h>
#include <src/pscf/pcompcabsfile.h>

using namespace std;
using namespace bagel;

//#define LOCAL_DEBUG

typedef shared_ptr<PMatrix1e> RefMatrix;
typedef shared_ptr<PGeometry> RefGeom;
typedef shared_ptr<PHcore> RefHcore;
typedef shared_ptr<PCoeff> RefCoeff;

#ifdef HAVE_LIBSLATER

pair<RefCoeff, RefCoeff> PMP2::generate_CABS() {

  // Form RI space which is a union of OBS and CABS.
  RefGeom newgeom(new PGeometry(*geom_));
  union_geom_ = newgeom;
  union_geom_->merge_obs_aux();

  shared_ptr<POverlap> union_overlap(new POverlap(union_geom_));
  shared_ptr<PTildeX> ri_coeff(new PTildeX(union_overlap));
  RefMatrix ri_reshaped(new PMatrix1e(coeff_, ri_coeff->ndim()));

  // SVD to project out OBS component. Note singular values are all 1 as OBS is a subset of RI space.
  RefMatrix uft(new PMatrix1e(union_overlap->ft()));
  uft->hermite();
  RefMatrix tmp(new PMatrix1e(*ri_coeff % *uft * *ri_reshaped));

  const int tmdim = tmp->mdim();
  const int tndim = tmp->ndim();

  shared_ptr<PMatrix1e> U, V;
  tie(U, V) = tmp->svd();

  RefMatrix Ured(new PMatrix1e(U, make_pair(tmdim, tndim)));
  RefCoeff coeff_cabs(new PCoeff(*ri_coeff * *Ured));
  coeff_cabs_ = coeff_cabs;

  RefMatrix coeff_fit(new PMatrix1e(coeff_, tndim));
  RefMatrix coeff_entire(new PMatrix1e(coeff_fit, coeff_cabs_));
  coeff_entire_ = coeff_entire;

  //(*coeff_entire_ % *uft * *coeff_entire_).rprint(15);

  pair<RefCoeff, RefCoeff> cabs_coeff_spl = coeff_cabs_->split(geom_->nbasis(), geom_->naux());

  return cabs_coeff_spl;
}


const std::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> PMP2::generate_hJ() {

  RefHcore uhc(new PHcore(union_geom_, true));
  RefMatrix coulombc = coulomb_runtime(true);

  RefMatrix ao_h(new PMatrix1e(uhc->ft()));
  RefMatrix ao_J(new PMatrix1e(coulombc->ft()));
  RefMatrix tmp(new PMatrix1e(*ao_h + *ao_J));
  ao_hJ_ = tmp;
  RefMatrix hJ(new PMatrix1e(*coeff_entire_ % *ao_hJ_ * *coeff_entire_));

  RefMatrix h_hJ_o(new PMatrix1e(hJ, make_pair(0, geom_->nbasis())));
  RefMatrix h_hJ_c(new PMatrix1e(hJ, make_pair(geom_->nbasis(), geom_->nbasis()+geom_->naux())));

  pair<RefMatrix, RefMatrix> h_hJ_o_pair = h_hJ_o->split(geom_->nbasis(), geom_->naux());
  pair<RefMatrix, RefMatrix> h_hJ_c_pair = h_hJ_c->split(geom_->nbasis(), geom_->naux());

  return make_tuple(h_hJ_o_pair.first, h_hJ_o_pair.second, h_hJ_c_pair.first, h_hJ_c_pair.second);
}


const std::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> PMP2::generate_K() {

  RefMatrix exchangec = exchange_runtime(true);
  RefMatrix ao_K(new PMatrix1e(exchangec->ft()));
  ao_K_ = ao_K;
  RefMatrix exchange(new PMatrix1e(*coeff_entire_ % *ao_K_ * *coeff_entire_));

  RefMatrix h_exchange_o(new PMatrix1e(exchange, make_pair(0, geom_->nbasis())));
  RefMatrix h_exchange_c(new PMatrix1e(exchange, make_pair(geom_->nbasis(), geom_->nbasis()+geom_->naux())));

  pair<RefMatrix, RefMatrix> h_exchange_o_pair = h_exchange_o->split(geom_->nbasis(), geom_->naux());
  pair<RefMatrix, RefMatrix> h_exchange_c_pair = h_exchange_c->split(geom_->nbasis(), geom_->naux());

  return make_tuple(h_exchange_o_pair.first, h_exchange_o_pair.second, h_exchange_c_pair.first, h_exchange_c_pair.second);
}


void PMP2::fill_in_cabs_matices() {
  // some preparation for P intermediate.

  // Hartree builder (needs modification!!!)
  // Index convention is *_upper_lower
  const tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> hJ4 = generate_hJ();
  // nbasis * nbasis size
  hJ_obs_obs_ = get<0>(hJ4);
  // nbasis * ncabs size
  hJ_obs_cabs_ = get<1>(hJ4);
  // ncabs * nbasis size
  hJ_cabs_obs_ = get<2>(hJ4);
  // ncabs * ncabs size
  hJ_cabs_cabs_ = get<3>(hJ4);
  // Exchange builder (needs modification!!!)
  const tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> K4 = generate_K();
  // nbasis * nbasis size
  K_obs_obs_ = get<0>(K4);
  // nbasis * ncabs size
  K_obs_cabs_ = get<1>(K4);
  // ncabs * nbasis size
  K_cabs_obs_ = get<2>(K4);
  // ncabs * ncabs size
  K_cabs_cabs_ = get<3>(K4);
  {
    RefMatrix fobs(new PMatrix1e(*hJ_obs_obs_ - *K_obs_obs_));
    RefMatrix fcabs(new PMatrix1e(*hJ_obs_cabs_ - *K_obs_cabs_));
    fock_obs_obs_ = fobs;
    fock_obs_cabs_ = fcabs;
  }
  {
    RefMatrix fobs(new PMatrix1e(*hJ_cabs_obs_ - *K_cabs_obs_));
    RefMatrix fcabs(new PMatrix1e(*hJ_cabs_cabs_ - *K_cabs_cabs_));
    fock_cabs_obs_ = fobs;
    fock_cabs_cabs_ = fcabs;
  }

}

#endif
