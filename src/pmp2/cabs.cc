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

using namespace std;
using namespace boost;

typedef shared_ptr<PMatrix1e> RefMatrix;
typedef shared_ptr<PGeometry> RefGeom;
typedef shared_ptr<PHcore> RefHcore;
typedef shared_ptr<PCoeff> RefCoeff;
typedef shared_ptr<PMOFile<std::complex<double> > > RefMOFile;

pair<RefCoeff, RefCoeff> PMP2::generate_CABS() {

  // Form RI space which is a union of OBS and CABS.
  RefGeom newgeom(new PGeometry(*geom_));
  union_geom_ = newgeom;
  union_geom_->merge_obs_cabs();

  shared_ptr<POverlap> union_overlap(new POverlap(union_geom_));
  shared_ptr<PTildeX> ri_coeff(new PTildeX(union_overlap));
  RefMatrix ri_reshaped(new PMatrix1e(coeff_, ri_coeff->ndim()));

  // SVD to project out OBS component. Note singular values are all 1 as OBS is a subset of RI space.
  RefMatrix tmp(new PMatrix1e(*ri_coeff % union_overlap->ft() * *ri_reshaped));

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

  pair<RefCoeff, RefCoeff> cabs_coeff_spl = coeff_cabs_->split(geom_->nbasis(), geom_->ncabs());

  return cabs_coeff_spl;
}


// Hartree-weighted (i.e., Fock except exchange) index space to be used in B intermediate evaluator.
RefMatrix PMP2::generate_hJ_obs(RefMOFile eri_Ip_Ip) {

  // Computes hcore in k-space.
  RefHcore hc(new PHcore(geom_));
  RefMatrix aohcore(new PMatrix1e(hc->ft()));
  RefMatrix mohcore(new PMatrix1e(*coeff_ % *aohcore * *coeff_));

  RefMatrix coulomb = eri_Ip_Ip->contract_density_J();
  RefMatrix hartree(new PMatrix1e(*mohcore + *coulomb));

  return hartree;
}


// Hartree-weighted (i.e., Fock except exchange) index space to be used in B intermediate evaluator.
RefMatrix PMP2::generate_hJ_cabs(RefMOFile eri_Ip_IA) {

  // Computes hcore in k-space.
  // Needs to pass that this PHcore is for union_geom_ (i.e. true in the second argument)
  // Assuming that the cost for this operation is negligible.
  RefHcore uhc(new PHcore(union_geom_, true));
  RefMatrix aohcore(new PMatrix1e(uhc->ft()));
  RefMatrix mohcore(new PMatrix1e(*coeff_entire_ % *aohcore * *coeff_entire_));

  RefMatrix mohcore_cabs = mohcore->split(geom_->nbasis(), geom_->ncabs()).second;
  RefMatrix mohcore_obs_cabs_block(new PMatrix1e(mohcore_cabs, make_pair(0, geom_->nbasis())));

  RefMatrix coulomb = eri_Ip_IA->contract_density_J();
  RefMatrix hartree(new PMatrix1e(*mohcore_obs_cabs_block + *coulomb));

  return hartree;
}

