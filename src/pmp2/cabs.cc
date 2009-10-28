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

typedef shared_ptr<PMatrix1e> RefMatrix;
typedef shared_ptr<PGeometry> RefGeom;
typedef shared_ptr<PHcore> RefHcore;
typedef shared_ptr<PCoeff> RefCoeff;
typedef shared_ptr<PMOFile<complex<double> > > RefMOFile;

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
RefMatrix PMP2::generate_hJ_obs_obs() {

  // TODO INEFFICIENT CODE!!! Hartree matrix needs to be constructed in AO basis.
  RefMOFile eri_pI_pI = eri_obs_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                               0, nbasis_, 0, nocc_,
                                               0, nbasis_, 0, nocc_, "h+J builder (OBS-OBS; pp)");
  eri_pI_pI->sort_inside_blocks();

  // Computes hcore in k-space.
  RefHcore hc(new PHcore(geom_));
  RefMatrix aohcore(new PMatrix1e(hc->ft()));
  RefMatrix mohcore(new PMatrix1e(*coeff_ % *aohcore * *coeff_));

  RefMatrix coulomb = eri_pI_pI->contract_density_J();
  RefMatrix hartree(new PMatrix1e(*mohcore + *coulomb));

  return hartree;
}


// Hartree-weighted (i.e., Fock except exchange) index space to be used in B intermediate evaluator.
RefMatrix PMP2::generate_hJ_obs_cabs() {

  // TODO INEFFICIENT CODE!!! Hartree matrix needs to be constructed in AO basis.
  RefMOFile eri_pI_pI = eri_obs_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                               0, nbasis_, 0, nocc_,
                                               0, ncabs_, 0, nocc_, "h+J builder (OBS-CABS; pp)");
  RefMOFile eri_pI_xI = eri_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                         0, nbasis_, 0, nocc_,
                                                         0, ncabs_, 0, nocc_, "h+J builder (OBS-CABS; px)");

  eri_pI_pI->sort_inside_blocks();
  eri_pI_xI->sort_inside_blocks();
  RefMOFile eri_pI_AI(new PMOFile<complex<double> >(*eri_pI_xI + *eri_pI_pI));

  // Computes hcore in k-space.
  // Needs to pass that this PHcore is for union_geom_ (i.e. true in the second argument)
  // Assuming that the cost for this operation is negligible.
  RefHcore uhc(new PHcore(union_geom_, true));
  RefMatrix aohcore(new PMatrix1e(uhc->ft()));
  RefMatrix mohcore(new PMatrix1e(*coeff_entire_ % *aohcore * *coeff_entire_));

  RefMatrix mohcore_cabs = mohcore->split(geom_->nbasis(), geom_->ncabs()).second;
  RefMatrix mohcore_obs_cabs_block(new PMatrix1e(mohcore_cabs, make_pair(0, geom_->nbasis())));

  RefMatrix coulomb = eri_pI_AI->contract_density_J();
  RefMatrix hartree(new PMatrix1e(*mohcore_obs_cabs_block + *coulomb));

  return hartree;
}


RefMatrix PMP2::generate_K_obs_obs() {
  // TODO INEFFICIENT CODE!!! Hartree matrix needs to be constructed in AO basis.
  RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                               0, nocc_, 0, nbasis_,
                                               0, nbasis_, 0, nocc_, "K builder (OBS-OBS; pp)");
  eri_Ip_pI->sort_inside_blocks();
  RefMatrix exchange = eri_Ip_pI->contract_density_K();
  return exchange;
}


RefMatrix PMP2::generate_K_obs_cabs() {
  // TODO INEFFICIENT CODE!!! Hartree matrix needs to be constructed in AO basis.
  RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                               0, nocc_, 0, nbasis_,
                                               0, ncabs_, 0, nocc_, "K builder (OBS-CABS; pp)");
  RefMOFile eri_Ip_xI = eri_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                         0, nocc_, 0, nbasis_,
                                                         0, ncabs_, 0, nocc_, "K builder (OBS-CABS; px)");
  eri_Ip_pI->sort_inside_blocks();
  eri_Ip_xI->sort_inside_blocks();
  RefMOFile eri_Ip_AI(new PMOFile<complex<double> >(*eri_Ip_xI + *eri_Ip_pI));
  RefMatrix exchange = eri_Ip_AI->contract_density_K();
  return exchange;
}


pair<RefMatrix, RefMatrix> PMP2::generate_K_cabs_pair() {
  const double gamma = geom_->gamma();
  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs_d(new PCompCABSFile<ERIBatch>(geom_, gamma, false, true, false, false, false, "ERI CABS(j)"));
  eri_cabs_d->store_integrals();
  eri_cabs_d->reopen_with_inout();

  // TODO INEFFICIENT CODE!!! Hartree matrix needs to be constructed in AO basis.
  // Note there is no bra-ket symmetry in periodic calculations.
  RefMatrix exchange1;
  {
    RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, cabs_obs_, coeff_, coeff_,
                                                 0, nocc_, 0, ncabs_,
                                                 0, nbasis_, 0, nocc_, "K builder (CABS-OBS; pp)");
    RefMOFile eri_Ix_pI = eri_cabs_d->mo_transform_cabs_aux(coeff_, cabs_aux_, coeff_, coeff_,
                                                            0, nocc_, 0, ncabs_,
                                                            0, nbasis_, 0, nocc_, "K builder (CABS-OBS; xp)");
    eri_Ip_pI->sort_inside_blocks();
    eri_Ix_pI->sort_inside_blocks();
    RefMOFile eri_IA_pI(new PMOFile<complex<double> >(*eri_Ix_pI + *eri_Ip_pI));
    exchange1 = eri_IA_pI->contract_density_K();
  }

  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs_t(new PCompCABSFile<ERIBatch>(geom_, gamma, false, true, true, false, false, "ERI CABS(ja)"));
  eri_cabs_t->store_integrals();
  eri_cabs_t->reopen_with_inout();

  // TODO INEFFICIENT CODE!!! Hartree matrix needs to be constructed in AO basis.
  // Note there is no bra-ket symmetry in periodic calculations.
  RefMatrix exchange2;
  {
    RefMOFile eri_Ix_xI = eri_cabs_t->mo_transform_cabs_aux(coeff_, cabs_aux_, cabs_aux_, coeff_,
                                                            0, nocc_, 0, ncabs_,
                                                            0, ncabs_, 0, nocc_, "K builder (CABS-CABS; xx)");
    RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, cabs_obs_, cabs_obs_, coeff_,
                                                 0, nocc_, 0, ncabs_,
                                                 0, ncabs_, 0, nocc_, "K builder (CABS-CABS; pp)");
    RefMOFile eri_Ix_pI = eri_cabs_d->mo_transform_cabs_aux(coeff_, cabs_aux_, cabs_obs_, coeff_,
                                                            0, nocc_, 0, ncabs_,
                                                            0, ncabs_, 0, nocc_, "K builder (CABS-OBS; xp)");
    RefMOFile eri_Ip_xI = eri_cabs_->mo_transform_cabs_aux(coeff_, cabs_obs_, cabs_aux_, coeff_,
                                                           0, nocc_, 0, ncabs_,
                                                           0, ncabs_, 0, nocc_, "K builder (OBS-CABS; px)");
    eri_Ip_pI->sort_inside_blocks();
    eri_Ix_pI->sort_inside_blocks();
    eri_Ip_xI->sort_inside_blocks();
    eri_Ix_xI->sort_inside_blocks();
    RefMOFile eri_IA_AI(new PMOFile<complex<double> >(*eri_Ip_pI + *eri_Ix_pI + *eri_Ip_xI + *eri_Ix_xI));
    exchange2 = eri_IA_AI->contract_density_K();
  }
  return make_pair(exchange1, exchange2);
}

