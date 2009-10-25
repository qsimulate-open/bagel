//
// Author : Toru Shiozaki
// Date   : August 2009
//

#include <cstring>
#include <iostream>
#include <algorithm>
#include <src/pmp2/pmp2.h>
#include <src/macros.h>
#include <src/scf/scf_macros.h>
#include <src/slater/slaterbatch.h>
#include <src/util/paircompfile.h>
#include <src/util/pcompcabsfile.h>
#include <src/util/pmofile.h>

typedef boost::shared_ptr<Atom> RefAtom;
typedef boost::shared_ptr<PGeometry> RefGeom;
typedef boost::shared_ptr<PCoeff> RefPCoeff;
typedef boost::shared_ptr<Shell> RefShell;
typedef boost::shared_ptr<PMOFile<std::complex<double> > > RefPMOFile;
typedef boost::shared_ptr<PMatrix1e> RefMatrix;

using namespace std;
using namespace boost;

PMP2::PMP2(const RefGeom g, const RefPCoeff co, const vector<double> eg, const shared_ptr<PCompFile<ERIBatch> > fl)
 : geom_(g), coeff_(co), eig_(eg.begin(), eg.end()), ao_eri_(fl) {

  nfrc_ = geom_->nfrc() / 2;
  nocc_ = geom_->nocc() / 2;
  nocc_act_ = nocc_ - nfrc_;
  nbasis_ = geom_->nbasis();
  nvir_ = nbasis_ - nocc_;
  noovv_ = nocc_act_ * nocc_act_ * nbasis_ * nbasis_;
  ncabs_ = geom_->ncabs();

  assert(geom_->ncabs() != 0);

}


PMP2::~PMP2() {

}


void PMP2::compute() {
  const double gamma = geom_->gamma();

  // AO ERI has been computed in the SCF class.

  cout << "  === Periodic MP2 calculation ===" << endl << endl;
  // Fully transform aa/ii integrals and dump them to disk (... focus is on MP2-R12).
  eri_ii_pp_ = ao_eri_->mo_transform(coeff_, nfrc_, nocc_, nfrc_, nocc_,
                                             0, nbasis_, 0, nbasis_, "ERI (pp/ii)");
  eri_ii_pp_->sort_inside_blocks();

  ///////////////////////////////////////////////
  // Compute the conventional MP2 contribution
  /////////////////////////////////////////////
  compute_conv_mp2();

  cout << "  === Periodic MP2-R12 calculation ===" << endl << endl;


  ////////// //////////////////////////////
  // Slater & Yukawa potential integrals
  ///////////////////////////////////////
  shared_ptr<PairCompFile<SlaterBatch> > stg_yp(new PairCompFile<SlaterBatch>(geom_, gamma, "Slater and Yukawa ints"));
  stg_yp->store_integrals();
  stg_yp->reopen_with_inout();
  shared_ptr<PCompFile<SlaterBatch> > stg = stg_yp->first();
  shared_ptr<PCompFile<SlaterBatch> > yp  = stg_yp->second();


  //////////////////////////////////////
  // Yukawa ii/ii for V intermediate
  ////////////////////////////////////
  RefPMOFile yp_ii_ii = yp->mo_transform(coeff_, nfrc_, nocc_, nfrc_, nocc_,
                                                 nfrc_, nocc_, nfrc_, nocc_, "Yukawa (ii/ii)");
  yp_ii_ii->sort_inside_blocks();
  RefPMOFile stg_ii_pp = stg->mo_transform(coeff_, nfrc_, nocc_, nfrc_, nocc_,
                                                   0, nbasis_, 0, nbasis_, "Slater (pp/ii)");
  stg_ii_pp->sort_inside_blocks();


  ////////////////////
  // CABS integrals
  //////////////////
  shared_ptr<PCompCABSFile<ERIBatch> >eri_cabs(new PCompCABSFile<ERIBatch>(geom_, gamma, false, "ERI CABS"));
  eri_cabs->store_integrals();
  eri_cabs->reopen_with_inout();

  shared_ptr<PCompCABSFile<SlaterBatch> >stg_cabs(new PCompCABSFile<SlaterBatch>(geom_, gamma, false, "Slater CABS"));
  stg_cabs->store_integrals();
  stg_cabs->reopen_with_inout();


  ///////////////////////////
  // Coefficients of CABS;
  /////////////////////////
  pair<RefMatrix, RefMatrix> cabs_pairs = generate_CABS();
  cabs_obs_ = cabs_pairs.first;
  cabs_aux_ = cabs_pairs.second;
  // Setting actual number of CABS.
  ncabs_ = cabs_obs_->mdim();


  ///////////////////////////////
  // MO transformation for ERI
  /////////////////////////////
  RefPMOFile eri_ii_ip = ao_eri_->mo_transform_cabs_obs(coeff_, cabs_obs_,
                                                        nfrc_, nocc_, nfrc_, nocc_,
                                                        0, nocc_, 0, ncabs_, "V^ia'_ii, OBS part");
  eri_ii_ip->sort_inside_blocks();

  RefPMOFile eri_ii_ix = eri_cabs->mo_transform_cabs_aux(coeff_, cabs_aux_,
                                                         nfrc_, nocc_, nfrc_, nocc_,
                                                         0, nocc_, 0, ncabs_, "V^ia'_ii, auxiliary functions");
  eri_ii_ix->sort_inside_blocks();
  RefPMOFile eri_ii_iA(new PMOFile<complex<double> >(*eri_ii_ix + *eri_ii_ip));

  ///////////////////////////////
  // MO transformation for STG
  /////////////////////////////
  RefPMOFile stg_ii_ip = stg->mo_transform_cabs_obs(coeff_, cabs_obs_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, nocc_, 0, ncabs_, "F^ia'_ii, OBS part");
  stg_ii_ip->sort_inside_blocks();
  RefPMOFile stg_ii_ix = stg_cabs->mo_transform_cabs_aux(coeff_, cabs_aux_,
                                                         nfrc_, nocc_, nfrc_, nocc_,
                                                         0, nocc_, 0, ncabs_, "F^ia'_ii, auxiliary functions");
  stg_ii_ix->sort_inside_blocks();
  RefPMOFile stg_ii_iA(new PMOFile<complex<double> >(*stg_ii_ix + *stg_ii_ip));


  //////////////////////
  // V intermediate
  ////////////////////
  RefPMOFile vF = stg_ii_pp->contract(eri_ii_pp_, nfrc_, nocc_, nfrc_, nocc_, "F * v (ii/ii) OBS");
  RefPMOFile V_obs(new PMOFile<complex<double> >(*yp_ii_ii - *vF));

  RefPMOFile V_cabs = stg_ii_iA->contract(eri_ii_iA, nfrc_, nocc_, nfrc_, nocc_, "F * v (ii/ii) CABS");

  RefPMOFile V(new PMOFile<complex<double> >(*V_obs - *V_cabs - *(V_cabs->flip())));
  complex<double> en_vt = V->get_energy_one_amp();


  ///////////////////////////////////////
  // Direct contribution to R12 energy
  /////////////////////////////////////
  cout << "  * Built V intermediate." << endl;
  cout << "  F12 energy (Vt): " << setprecision(10) << en_vt.real() << endl << endl;


  /////////////////////////////////////
  // Slater & Yukawa ints with 2gamma
  ///////////////////////////////////
  shared_ptr<PairCompFile<SlaterBatch> >
    stg_yp2(new PairCompFile<SlaterBatch>(geom_, 2.0 * gamma, "Slater and Yukawa ints (2gamma)"));
  stg_yp2->store_integrals();
  stg_yp2->reopen_with_inout();
  shared_ptr<PCompFile<SlaterBatch> > stg2 = stg_yp2->first();
  shared_ptr<PCompFile<SlaterBatch> > yp2  = stg_yp2->second();


  ////////////////////
  // X intermediate
  //////////////////
  RefPMOFile stg2_ii_ii = stg2->mo_transform(coeff_, nfrc_, nocc_, nfrc_, nocc_,
                                                     nfrc_, nocc_, nfrc_, nocc_, "Slater (ii/ii) 2gamma");
  stg2_ii_ii->sort_inside_blocks();
  RefPMOFile FF = stg_ii_pp->contract(stg_ii_pp, nfrc_, nocc_, nfrc_, nocc_, "F * F (ii/ii) OBS");
  RefPMOFile X_obs(new PMOFile<complex<double> >(*stg2_ii_ii - *FF));

  RefPMOFile X_cabs = stg_ii_iA->contract(stg_ii_iA, nfrc_, nocc_, nfrc_, nocc_, "F * F (ii/ii) CABS");
  RefPMOFile X(new PMOFile<complex<double> >(*X_obs - *X_cabs - *(X_cabs->flip())));


  /////////////////////
  // B intermediate
  ///////////////////

  // T intermediate (direct)
  RefPMOFile T(new PMOFile<complex<double> >(*stg2_ii_ii * (gamma*gamma)));
  T->rprint();

  // Q intermediate (made of X * h)
  pair<RefMatrix, RefMatrix> hj = generate_hJ();

#ifdef LOCAL_DEBUG_PMP2
  V->print();
  cout << endl;
  X->print();
#endif

}


