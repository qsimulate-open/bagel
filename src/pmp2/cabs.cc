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
  RefMatrix ri_reshaped(new PMatrix1e(coeff_, ri_coeff->ndim(), coeff_->mdim()));

  // SVD to project out OBS component. Note singular values are all 1 as OBS is a subset of RI space.
  RefMatrix tmp(new PMatrix1e(*ri_coeff % union_overlap->ft() * *ri_reshaped));
  RefMatrix U(new PMatrix1e(geom_, tmp->ndim(), tmp->ndim()));
  RefMatrix V(new PMatrix1e(geom_, tmp->mdim(), tmp->mdim()));
  tmp->svd(U, V);

  RefMatrix Ured(new PMatrix1e(U, make_pair(tmp->mdim(), tmp->ndim())));
  RefCoeff cabs_coeff(new PCoeff(*ri_coeff * *Ured));

  pair<RefCoeff, RefCoeff> cabs_coeff_spl = cabs_coeff->split(geom_->nbasis(), geom_->ncabs());

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
pair<RefMatrix, RefMatrix> PMP2::generate_hJ_iA(RefMOFile eri_Ii_IA) {

  // Computes hcore in k-space.
  RefHcore union_hcore(new PHcore(union_geom_));
  RefMatrix hcore(new PMatrix1e(union_hcore->ft()));

//hcore->rprint(4);

  RefMatrix coulomb = eri_Ii_IA->contract_density_J();

  coulomb->rprint(4);

  pair<RefMatrix, RefMatrix> out;
  return out;
}

