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


const boost::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> PMP2::generate_hJ() const {

#define USE_BUILDER
  // hcore contribution.
  RefHcore uhc(new PHcore(union_geom_, true));

#ifdef USE_BUILDER

  RefMatrix ao_hJ(new PMatrix1e((*uhc + *coulomb_runtime()).ft()));
  RefMatrix hJ(new PMatrix1e(*coeff_entire_ % *ao_hJ * *coeff_entire_));

#else
  RefMatrix aohcore(new PMatrix1e(uhc->ft()));
  RefMatrix mohcore(new PMatrix1e(*coeff_entire_ % *aohcore * *coeff_entire_));
  // coulomb obs_obs block
  RefMatrix coulomb_oo;
  {
    RefMOFile eri_pI_pI = eri_obs_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                                 0, nbasis_, 0, nocc_,
                                                 0, nbasis_, 0, nocc_, "h+J builder (OBS-OBS; pp)");


    coulomb_oo = eri_pI_pI->contract_density_J();
  }

  // coulomb obs_cabs block
  RefMatrix coulomb_oc;
  {
    RefMOFile eri_pI_pI = eri_obs_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                                 0, nbasis_, 0, nocc_,
                                                 0, ncabs_, 0, nocc_, "h+J builder (OBS-CABS; pp)");
    RefMOFile eri_pI_xI = eri_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                           0, nbasis_, 0, nocc_,
                                                           0, ncabs_, 0, nocc_, "h+J builder (OBS-CABS; px)");
    RefMOFile eri_pI_AI(new PMOFile<complex<double> >(*eri_pI_xI + *eri_pI_pI));

    coulomb_oc = eri_pI_AI->contract_density_J();
  }


  const double gamma = geom_->gamma();
  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs_d(new PCompCABSFile<ERIBatch>(geom_, gamma, true, false, false, false, false, "ERI CABS(j)"));
  eri_cabs_d->store_integrals();
  eri_cabs_d->reopen_with_inout();

  // coulomb cabs_obs block
  RefMatrix coulomb_co;
  {
    RefMOFile eri_AI_pI = eri_obs_->mo_transform(cabs_obs_, coeff_, coeff_, coeff_,
                                                 0, ncabs_, 0, nocc_,
                                                 0, nbasis_, 0, nocc_, "hJ builder (CABS-OBS; pp)");
    *eri_AI_pI += *(eri_cabs_d->mo_transform_cabs_aux(cabs_aux_, coeff_, coeff_, coeff_,
                                                      0, ncabs_, 0, nocc_,
                                                      0, nbasis_, 0, nocc_, "hJ builder (CABS-OBS; xp)"));
    coulomb_co = eri_AI_pI->contract_density_J();
  }

  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs_t(new PCompCABSFile<ERIBatch>(geom_, gamma, true, false, true, false, false, "ERI CABS(ja)"));
  eri_cabs_t->store_integrals();
  eri_cabs_t->reopen_with_inout();

  // coulomb cabs_cabs block
  RefMatrix coulomb_cc;
  {
    RefMOFile eri_AI_AI = eri_cabs_t->mo_transform_cabs_aux(cabs_aux_, coeff_, cabs_aux_, coeff_,
                                                            0, ncabs_, 0, nocc_,
                                                            0, ncabs_, 0, nocc_, "hJ builder (CABS-CABS; xx)");
    *eri_AI_AI += *(eri_obs_->mo_transform(cabs_obs_, coeff_, cabs_obs_, coeff_,
                                           0, ncabs_, 0, nocc_,
                                           0, ncabs_, 0, nocc_, "hJ builder (CABS-CABS; pp)"));
    *eri_AI_AI += *(eri_cabs_d->mo_transform_cabs_aux(cabs_aux_, coeff_, cabs_obs_, coeff_,
                                                      0, ncabs_, 0, nocc_,
                                                      0, ncabs_, 0, nocc_, "hJ builder (CABS-CABS; xp)"));
    *eri_AI_AI += *(eri_cabs_->mo_transform_cabs_aux(cabs_obs_, coeff_, cabs_aux_, coeff_,
                                                     0, ncabs_, 0, nocc_,
                                                     0, ncabs_, 0, nocc_, "hJ builder (CABS-CABS; px)"));
    coulomb_cc = eri_AI_AI->contract_density_J();
  }

  RefMatrix coulomb(new PMatrix1e(coulomb_oo->merge(coulomb_oc), coulomb_co->merge(coulomb_cc)));
  coulomb->conj();
  RefMatrix hJ(new PMatrix1e(*coulomb+*mohcore));
#endif

  RefMatrix h_hJ_o(new PMatrix1e(hJ, make_pair(0, geom_->nbasis())));
  RefMatrix h_hJ_c(new PMatrix1e(hJ, make_pair(geom_->nbasis(), geom_->nbasis()+geom_->ncabs())));

  pair<RefMatrix, RefMatrix> h_hJ_o_pair = h_hJ_o->split(geom_->nbasis(), geom_->ncabs());
  pair<RefMatrix, RefMatrix> h_hJ_c_pair = h_hJ_c->split(geom_->nbasis(), geom_->ncabs());

  return make_tuple(h_hJ_o_pair.first, h_hJ_o_pair.second, h_hJ_c_pair.first, h_hJ_c_pair.second);
}




const boost::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> PMP2::generate_K() const {
  // TODO INEFFICIENT CODE!!! Hartree matrix needs to be constructed in AO basis.

  // obs_obs block
  RefMatrix exchange_oo;
  {
    RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                                 0, nocc_, 0, nbasis_,
                                                 0, nbasis_, 0, nocc_, "K builder (OBS-OBS; pp)");
    exchange_oo = eri_Ip_pI->contract_density_K();
  }

  // obs_cabs block
  RefMatrix exchange_oc;
  {
    RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                                 0, nocc_, 0, nbasis_,
                                                 0, ncabs_, 0, nocc_, "K builder (OBS-CABS; pp)");
    RefMOFile eri_Ip_xI = eri_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                           0, nocc_, 0, nbasis_,
                                                           0, ncabs_, 0, nocc_, "K builder (OBS-CABS; px)");
    RefMOFile eri_Ip_AI(new PMOFile<complex<double> >(*eri_Ip_xI + *eri_Ip_pI));
    exchange_oc = eri_Ip_AI->contract_density_K();
  }

  const double gamma = geom_->gamma();
  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs_d(new PCompCABSFile<ERIBatch>(geom_, gamma, false, true, false, false, false, "ERI CABS(j)"));
  eri_cabs_d->store_integrals();
  eri_cabs_d->reopen_with_inout();

  //  cabs_obs block
  RefMatrix exchange_co;
  {
    RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, cabs_obs_, coeff_, coeff_,
                                                 0, nocc_, 0, ncabs_,
                                                 0, nbasis_, 0, nocc_, "K builder (CABS-OBS; pp)");
    RefMOFile eri_Ix_pI = eri_cabs_d->mo_transform_cabs_aux(coeff_, cabs_aux_, coeff_, coeff_,
                                                            0, nocc_, 0, ncabs_,
                                                            0, nbasis_, 0, nocc_, "K builder (CABS-OBS; xp)");
    RefMOFile eri_IA_pI(new PMOFile<complex<double> >(*eri_Ix_pI + *eri_Ip_pI));
    exchange_co = eri_IA_pI->contract_density_K();
  }

  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs_t(new PCompCABSFile<ERIBatch>(geom_, gamma, false, true, true, false, false, "ERI CABS(ja)"));
  eri_cabs_t->store_integrals();
  eri_cabs_t->reopen_with_inout();

  // cabs-cabs block
  RefMatrix exchange_cc;
  {
    RefMOFile eri_Ix_xI = eri_cabs_t->mo_transform_cabs_aux(coeff_, cabs_aux_, cabs_aux_, coeff_,
                                                            0, nocc_, 0, ncabs_,
                                                            0, ncabs_, 0, nocc_, "K builder (CABS-CABS; xx)");
    RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, cabs_obs_, cabs_obs_, coeff_,
                                                 0, nocc_, 0, ncabs_,
                                                 0, ncabs_, 0, nocc_, "K builder (CABS-CABS; pp)");
    RefMOFile eri_Ix_pI = eri_cabs_d->mo_transform_cabs_aux(coeff_, cabs_aux_, cabs_obs_, coeff_,
                                                            0, nocc_, 0, ncabs_,
                                                            0, ncabs_, 0, nocc_, "K builder (CABS-CABS; xp)");
    RefMOFile eri_Ip_xI = eri_cabs_->mo_transform_cabs_aux(coeff_, cabs_obs_, cabs_aux_, coeff_,
                                                           0, nocc_, 0, ncabs_,
                                                           0, ncabs_, 0, nocc_, "K builder (CABS-CABS; px)");
    RefMOFile eri_IA_AI(new PMOFile<complex<double> >(*eri_Ip_pI + *eri_Ix_pI + *eri_Ip_xI + *eri_Ix_xI));
    exchange_cc = eri_IA_AI->contract_density_K();
  }

  // first form entire matrix
  RefMatrix exchange(new PMatrix1e(exchange_oo->merge(exchange_oc), exchange_co->merge(exchange_cc)));
  // this is an important step!!!
  exchange->conj();

  RefMatrix h_exchange_o(new PMatrix1e(exchange, make_pair(0, geom_->nbasis())));
  RefMatrix h_exchange_c(new PMatrix1e(exchange, make_pair(geom_->nbasis(), geom_->nbasis()+geom_->ncabs())));

  pair<RefMatrix, RefMatrix> h_exchange_o_pair = h_exchange_o->split(geom_->nbasis(), geom_->ncabs());
  pair<RefMatrix, RefMatrix> h_exchange_c_pair = h_exchange_c->split(geom_->nbasis(), geom_->ncabs());

  return make_tuple(h_exchange_o_pair.first, h_exchange_o_pair.second, h_exchange_c_pair.first, h_exchange_c_pair.second);
}


