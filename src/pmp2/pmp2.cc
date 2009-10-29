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
#include <src/util/paircompfile.h>
#include <src/util/pcompcabsfile.h>
#include <src/util/pmofile.h>

typedef boost::shared_ptr<Atom> RefAtom;
typedef boost::shared_ptr<PGeometry> RefGeom;
typedef boost::shared_ptr<PCoeff> RefPCoeff;
typedef boost::shared_ptr<Shell> RefShell;
typedef boost::shared_ptr<PMOFile<std::complex<double> > > RefPMOFile;
typedef boost::shared_ptr<PMatrix1e> RefMatrix;

// TODO I have not symmetrize intermediates to Hermitian as we are now using fixed amplitudes.

using namespace std;
using namespace boost;

PMP2::PMP2(const RefGeom g, const RefPCoeff co, const vector<double> eg, const shared_ptr<PCompFile<ERIBatch> > fl)
 : geom_(g), coeff_(co), eig_(eg.begin(), eg.end()), eri_obs_(fl) {

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
  eri_ii_pp_ = eri_obs_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                      nfrc_, nocc_, nfrc_, nocc_,
                                      0, nbasis_, 0, nbasis_, "ERI (pp/ii)");

  ///////////////////////////////////////////////
  // Compute the conventional MP2 contribution
  /////////////////////////////////////////////
  compute_conv_mp2();

  cout << "  === Periodic MP2-R12 calculation ===" << endl << endl;


  ////////// //////////////////////////////
  // Slater & Yukawa potential integrals
  ///////////////////////////////////////
  {
    shared_ptr<PairCompFile<SlaterBatch> > stg_yp(new PairCompFile<SlaterBatch>(geom_, gamma, "Slater and Yukawa ints"));
    stg_yp->store_integrals();
    stg_yp->reopen_with_inout();
    shared_ptr<PCompFile<SlaterBatch> > stg = stg_yp->first();
    shared_ptr<PCompFile<SlaterBatch> > yp  = stg_yp->second();
    stg_ = stg;
    yp_ = yp;
  }


  //////////////////////////////////////
  // Yukawa ii/ii for V intermediate
  ////////////////////////////////////
  {
    RefPMOFile yp_ii_ii = yp_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            nfrc_, nocc_, nfrc_, nocc_, "Yukawa (ii/ii)");
    yp_ii_ii_ = yp_ii_ii;
    RefPMOFile stg_ii_pp = stg_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              0, nbasis_, 0, nbasis_, "Slater (pp/ii)");
    stg_ii_pp_ = stg_ii_pp;
  }


  ////////////////////
  // CABS integrals
  //////////////////
  {
    shared_ptr<PCompCABSFile<ERIBatch> >
      eri_cabs(new PCompCABSFile<ERIBatch>(geom_, gamma, false, false, true, false, false, "ERI CABS"));
    eri_cabs->store_integrals();
    eri_cabs->reopen_with_inout();
    eri_cabs_ = eri_cabs;
  }
  {
    shared_ptr<PCompCABSFile<SlaterBatch> >
      stg_cabs(new PCompCABSFile<SlaterBatch>(geom_, gamma, false, false, true, false, false, "Slater CABS"));
    stg_cabs->store_integrals();
    stg_cabs->reopen_with_inout();
    stg_cabs_ = stg_cabs;
  }

  ///////////////////////////
  // Coefficients of CABS;
  /////////////////////////
  pair<RefCoeff, RefCoeff> cabs_pairs = generate_CABS();
  cabs_obs_ = cabs_pairs.first;
  cabs_aux_ = cabs_pairs.second;
  // Setting actual number of CABS.
  ncabs_ = cabs_obs_->mdim();


  ///////////////////////////////
  // MO transformation for ERI
  /////////////////////////////
  {
    RefPMOFile eri_ii_pi = eri_obs_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                                  nfrc_, nocc_, nfrc_, nocc_,
                                                  0, ncabs_, 0, nocc_, "v^ia'_ii, OBS part");

    RefPMOFile eri_ii_xi = eri_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                            nfrc_, nocc_, nfrc_, nocc_,
                                                            0, ncabs_, 0, nocc_, "v^ia'_ii, auxiliary functions");
    RefPMOFile eri_ii_Ai(new PMOFile<complex<double> >(*eri_ii_xi + *eri_ii_pi));
    eri_ii_Ai_ = eri_ii_Ai;
  }

  ///////////////////////////////
  // MO transformation for STG
  /////////////////////////////
  {
    RefPMOFile stg_ii_pi = stg_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              0, ncabs_, 0, nocc_, "F^ia'_ii, OBS part");
    RefPMOFile stg_ii_xi = stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                            nfrc_, nocc_, nfrc_, nocc_,
                                                            0, ncabs_, 0, nocc_, "F^ia'_ii, auxiliary functions");
    RefPMOFile stg_ii_Ai(new PMOFile<complex<double> >(*stg_ii_xi + *stg_ii_pi));
    stg_ii_Ai_ = stg_ii_Ai;
  }


  //////////////////////
  // V intermediate
  ////////////////////
  {
    RefPMOFile vF = stg_ii_pp_->contract(eri_ii_pp_, "F * v (ii/ii) OBS");
    RefPMOFile V_obs(new PMOFile<complex<double> >(*yp_ii_ii_ - *vF));

    RefPMOFile V_cabs = stg_ii_Ai_->contract(eri_ii_Ai_, "F * v (ii/ii) CABS");

    V_cabs->flip_symmetry();
    RefPMOFile V_pre(new PMOFile<complex<double> >(*V_obs - *V_cabs));

    V_ = V_pre;
  }
  complex<double> en_vt = V_->get_energy_one_amp();


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

  shared_ptr<PCompCABSFile<SlaterBatch> >
    stg2_cabs(new PCompCABSFile<SlaterBatch>(geom_, 2.0 * gamma, false, false, true, false,
                                             false, "Slater CABS (2gamma)"));
  stg2_cabs->store_integrals();
  stg2_cabs->reopen_with_inout();


  ////////////////////
  // X intermediate
  //////////////////
  {
    RefPMOFile stg2_ii_ii = stg2->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                               nfrc_, nocc_, nfrc_, nocc_,
                                               nfrc_, nocc_, nfrc_, nocc_, "Slater (ii/ii) 2gamma");
    stg2_ii_ii_ = stg2_ii_ii;
    RefPMOFile FF = stg_ii_pp_->contract(stg_ii_pp_, "F * F (ii/ii) OBS");
    RefPMOFile X_obs(new PMOFile<complex<double> >(*stg2_ii_ii - *FF));

    RefPMOFile X_cabs = stg_ii_Ai_->contract(stg_ii_Ai_, "F * F (ii/ii) CABS");

    X_cabs->flip_symmetry();
    RefPMOFile X_pre(new PMOFile<complex<double> >(*X_obs - *X_cabs));

    X_ = X_pre;
  }


  /////////////////////
  // B intermediate
  ///////////////////
  {
    // T intermediate (direct)
    RefPMOFile T(new PMOFile<complex<double> >(*stg2_ii_ii_ * (gamma*gamma)));

    // Q intermediate (made of X * h)
    RefPMOFile Q;
    {
      // Hartree builder (needs modification!!!)
      // nbasis * nbasis size
      hJ_obs_obs_ = generate_hJ_obs_obs();
      // nbasis * ncabs size
      hJ_obs_cabs_ = generate_hJ_obs_cabs();


      // Hartree weighted index
      RefMatrix hj_ip(new PMatrix1e(hJ_obs_obs_, make_pair(0, nocc_)));
      RefPCoeff chj_ip(new PCoeff(*coeff_ * *hj_ip));

      RefMatrix hj_iA(new PMatrix1e(hJ_obs_cabs_, make_pair(0, nocc_)));
      RefPCoeff chj_iA_comb(new PCoeff(*coeff_cabs_ * *hj_iA));

      pair<RefMatrix, RefMatrix> chj_iA = chj_iA_comb->split(geom_->nbasis(), geom_->ncabs());
      RefMatrix chj_iA_obs = chj_iA.first;
      RefCoeff chj_iA_cabs(new PCoeff(*chj_iA.second));

      *chj_ip += *chj_iA_obs;
      chj_ip->scale(0.5);
      chj_iA_cabs->scale(0.5);

      // MO transform using Hartree-weighted index
      RefPMOFile X_ii_ih = stg2->mo_transform(coeff_, coeff_,  chj_ip, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              nfrc_, nocc_, nfrc_, nocc_, "Q intermediate: stg2 (OBS) 1/2");
      RefPMOFile X_ii_ih_cabs = stg2_cabs->mo_transform_cabs_aux(coeff_, coeff_,  chj_iA_cabs, coeff_,
                                                                 nfrc_, nocc_, nfrc_, nocc_,
                                                                 nfrc_, nocc_, nfrc_, nocc_, "Q intermediate: stg2 (CABS) 1/2");

      *X_ii_ih += *X_ii_ih_cabs;

      X_ii_ih->flip_symmetry();
      RefPMOFile Qtmp(new PMOFile<complex<double> >(*X_ii_ih));
      Q = Qtmp;
      Q->scale(2.0);
    } // end of Q intermediate construction.

    Q->rprint();

    // some preparation for P intermediate.
    {
      // Exchange builder (needs modification!!!)
      // nbasis * nbasis size
      K_obs_obs_ = generate_K_obs_obs();
      // nbasis * ncabs size
      K_obs_cabs_ = generate_K_obs_cabs();

      pair<RefMatrix, RefMatrix> p = generate_K_cabs_pair();
      // ncabs * nbasis size
      K_cabs_obs_ = p.first;
      // ncabs * ncabs size
      K_cabs_cabs_ = p.second;

      RefMatrix fobs(new PMatrix1e(*hJ_obs_obs_ - *K_obs_obs_));
      fock_obs_obs_ = fobs;
    }
    // P intermediate R^PQ_ij K^R_P R^kl_RQ
    {
#if 0
      // construct entire matrix:
      RefMatrix K_obs_ri(new PMatrix1e(K_obs_obs_, K_cabs_obs_));
      RefMatrix K_cabs_ri(new PMatrix1e(K_obs_cabs_, K_cabs_cabs_));
      RefMatrix K_entire = K_obs_ri->merge(K_cabs_ri);
      coeff_entire_->rprint();
      RefPCoeff cK_AA(new PCoeff(*coeff_entire_ * *K_entire));
      cK_AA->rprint();

      shared_ptr<PCompCABSFile<SlaterBatch> >
        slater_cabs2(new PCompCABSFile<SlaterBatch>(geom_, gamma, false, false, true, true, false, "Slater CABS(ab)"));
      slater_cabs2->store_integrals();
      slater_cabs2->reopen_with_inout();
      RefMOFile slater_ii_xx = slater_cabs2->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cabs_aux_,
                                                                   nfrc_, nocc_, nfrc_, nocc_,
                                                                   0, ncabs_, 0, ncabs_, "Slater (CABS ab); xx)");
      RefMOFile slater_ii_xp = stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cabs_obs_,
                                                                nfrc_, nocc_, nfrc_, nocc_,
                                                                0, ncabs_, 0, ncabs_, "Slater (CABS ab); xp)");
      sort!
      slater_ii_xp->flip_symmetry();
      *slater_ii_xx += *slater_ii_xp;
      RefMOFile slater_ii_ = stg_->mo_transform(coeff_, coeff_, cabs_obs_, cabs_obs_,
                                                  nfrc_, nocc_, nfrc_, nocc_,
                                                  0, ncabs_, 0, ncabs_, "Slater (CABS ab); xx)");
#endif
    }

    RefCoeff cfobs(new PCoeff(*coeff_ * *fock_obs_obs_));
    // P3 intermediate R^m
    {
       RefMOFile p3_1 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, cfobs,
                                           nfrc_, nocc_, nfrc_, nocc_,
                                           0, ncabs_, 0, nocc_, "P3: R^Am_ij f^m_n R^kl_An, 1/2, OBS");
       *p3_1 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cfobs,
                                                   nfrc_, nocc_, nfrc_, nocc_,
                                                   0, ncabs_, 0, nocc_, "P3: R^mA_ij f^m_n R^kl_nA, 1/2, CABS"));
       RefMOFile p3_2 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                           nfrc_, nocc_, nfrc_, nocc_,
                                           0, ncabs_, 0, nocc_, "P3: R^Am_ij f^m_n R^kl_An, 1/2, OBS");
       *p3_2 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                   nfrc_, nocc_, nfrc_, nocc_,
                                                   0, ncabs_, 0, nocc_, "P3: R^mA_ij f^m_n R^kl_nA, 1/2, CABS"));
       RefPMOFile p3 = p3_1->contract(p3_2, "P3: R^Am_ij f^m_n R^kl_An");
       p3->flip_symmetry();
       p3->rprint();
     }

    // P4 intermediate R^pb_ij f^r_p R^ij_rb
    {
      RefMOFile p4_1 = stg_->mo_transform(coeff_, coeff_, cfobs, coeff_,
                                          nfrc_, nocc_, nfrc_, nocc_,
                                          0, nbasis_, nocc_ , nbasis_, "P4: R^pb_ij f^r_p R^kl_rb, 1/2");
      RefMOFile p4_2 = stg_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                          nfrc_, nocc_, nfrc_, nocc_,
                                          0, nbasis_, nocc_ , nbasis_, "P4: R^pb_ij f^r_p R^kl_rb, 2/2");
      RefPMOFile p4 = p4_1->contract(p4_2, "P4: R^pb_ij f^r_p R^ij_rb");
      p4->flip_symmetry();
      p4->rprint();

    }

  } // end of B intermediate construction.

#ifdef LOCAL_DEBUG_PMP2
  V->print();
  cout << endl;
  X->print();
#endif

}


