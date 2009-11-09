/*
 * cabs.cc
 *
 *  Created on: Oct 22, 2009
 *      Author: shiozaki
 */

#include <src/pmp2/pmp2.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/ptildex.h>
#include <src/pscf/phcore.h>
#include <src/pscf/pcoeff.h>
#include <src/util/pcompcabsfile.h>

using namespace std;
using namespace boost;

//#define LOCAL_DEBUG

typedef shared_ptr<PMatrix1e> RefMatrix;
typedef shared_ptr<PGeometry> RefGeom;
typedef shared_ptr<PHcore> RefHcore;
typedef shared_ptr<PCoeff> RefCoeff;

pair<RefCoeff, RefCoeff> PMP2::generate_CABS() {

  // Form RI space which is a union of OBS and CABS.
  RefGeom newgeom(new PGeometry(*geom_));
  union_geom_ = newgeom;
  union_geom_->merge_obs_cabs();

  shared_ptr<POverlap> union_overlap(new POverlap(union_geom_));
  shared_ptr<PTildeX> ri_coeff(new PTildeX(union_overlap));
  RefMatrix ri_reshaped(new PMatrix1e(coeff_, ri_coeff->ndim()));

  // SVD to project out OBS component. Note singular values are all 1 as OBS is a subset of RI space.
  RefMatrix uft(new PMatrix1e(union_overlap->ft()));
  uft->hermite();
  RefMatrix tmp(new PMatrix1e(*ri_coeff % *uft * *ri_reshaped));

  const int tmdim = tmp->mdim();
  const int tndim = tmp->ndim();

  RefMatrix U(new PMatrix1e(geom_, tndim, tndim));
  RefMatrix V(new PMatrix1e(geom_, tmdim, tmdim));
  tmp->svd(U, V);

  RefMatrix Ured(new PMatrix1e(U, make_pair(tmdim, tndim)));
  RefCoeff coeff_cabs(new PCoeff(*ri_coeff * *Ured));
  coeff_cabs_ = coeff_cabs;

  RefMatrix coeff_fit(new PMatrix1e(coeff_, tndim));
  RefMatrix coeff_entire(new PMatrix1e(coeff_fit, coeff_cabs_));
  coeff_entire_ = coeff_entire;

  //(*coeff_entire_ % *uft * *coeff_entire_).rprint(15);

  pair<RefCoeff, RefCoeff> cabs_coeff_spl = coeff_cabs_->split(geom_->nbasis(), geom_->ncabs());

  return cabs_coeff_spl;
}


const boost::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> PMP2::generate_hJ() {

  RefHcore uhc(new PHcore(union_geom_, true));
  RefMatrix coulombc = coulomb_runtime(true);

  RefMatrix ao_h(new PMatrix1e(uhc->ft()));
  RefMatrix ao_J(new PMatrix1e(coulombc->ft()));
  RefMatrix tmp(new PMatrix1e(*ao_h + *ao_J));
  ao_hJ_ = tmp;
  RefMatrix hJ(new PMatrix1e(*coeff_entire_ % *ao_hJ_ * *coeff_entire_));

  RefMatrix h_hJ_o(new PMatrix1e(hJ, make_pair(0, geom_->nbasis())));
  RefMatrix h_hJ_c(new PMatrix1e(hJ, make_pair(geom_->nbasis(), geom_->nbasis()+geom_->ncabs())));

  pair<RefMatrix, RefMatrix> h_hJ_o_pair = h_hJ_o->split(geom_->nbasis(), geom_->ncabs());
  pair<RefMatrix, RefMatrix> h_hJ_c_pair = h_hJ_c->split(geom_->nbasis(), geom_->ncabs());

  return make_tuple(h_hJ_o_pair.first, h_hJ_o_pair.second, h_hJ_c_pair.first, h_hJ_c_pair.second);
}


const boost::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> PMP2::generate_K() {

  RefMatrix exchangec = exchange_runtime(true);
  RefMatrix ao_K(new PMatrix1e(exchangec->ft()));
  ao_K_ = ao_K;
  RefMatrix exchange(new PMatrix1e(*coeff_entire_ % *ao_K_ * *coeff_entire_));

  RefMatrix h_exchange_o(new PMatrix1e(exchange, make_pair(0, geom_->nbasis())));
  RefMatrix h_exchange_c(new PMatrix1e(exchange, make_pair(geom_->nbasis(), geom_->nbasis()+geom_->ncabs())));

  pair<RefMatrix, RefMatrix> h_exchange_o_pair = h_exchange_o->split(geom_->nbasis(), geom_->ncabs());
  pair<RefMatrix, RefMatrix> h_exchange_c_pair = h_exchange_c->split(geom_->nbasis(), geom_->ncabs());

  return make_tuple(h_exchange_o_pair.first, h_exchange_o_pair.second, h_exchange_c_pair.first, h_exchange_c_pair.second);
}
