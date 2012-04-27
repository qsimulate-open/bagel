//
// Newint - Parallel electron correlation program.
// Filename: pmp2.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <cstring>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <tuple>
#include <src/pmp2/pmp2.h>
#include <src/macros.h>
#include <src/util/paircompfile.h>
#include <src/util/pcompcabsfile.h>
#include <src/util/pmofile.h>

typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<PGeometry> RefGeom;
typedef std::shared_ptr<PCoeff> RefPCoeff;
typedef std::shared_ptr<Shell> RefShell;
typedef std::shared_ptr<PMatrix1e> RefMatrix;

extern "C" { void start_up_slater_(); };

// TODO I have not symmetrize intermediates to Hermitian as we are now using fixed amplitudes.

#define DEBUG_PRINT

//#define ONLY_V_X
//#define ONLY_X

//#define ONLY_B
//#define ONLY_T_Q

//#define ONLY_P
//#define ONLY_P1
//#define SKIP_CASE1
//#define ONLY_P2
//#define ONLY_P3

using namespace std;

PMP2::PMP2(const RefGeom g, const RefPCoeff co, const vector<double> eg, const shared_ptr<PCompFile<ERIBatch> > fl, const bool hy2)
 : geom_(g), coeff_(co), eig_(eg.begin(), eg.end()), eri_obs_(fl), use_hy2_(hy2) {

  cout << "  === Periodic MP2 calculation ===" << endl << endl;

  cout << "    >> Gamma = " << geom_->gamma() << " is used in the MP2-R12 calculation." << endl << endl;
  if (use_hy2_) {
    cout << "    >> Hybrid approximation HY2 is invoked." << endl << endl;
  }

  nfrc_ = geom_->nfrc() / 2;
  nocc_ = geom_->nocc() / 2;
  nocc_act_ = nocc_ - nfrc_;
  nbasis_ = geom_->nbasis();
  nvir_ = nbasis_ - nocc_;
  noovv_ = nocc_act_ * nocc_act_ * nbasis_ * nbasis_;
  ncabs_ = geom_->naux();

  if (geom_->naux() == 0) {
    throw runtime_error("CABS should be specified.");
  }

  start_up_slater_();

}


PMP2::~PMP2() {

}


void PMP2::compute() {
  const double gamma = geom_->gamma();

  ///////////////////////////
  // Coefficients of CABS;
  /////////////////////////
  pair<RefCoeff, RefCoeff> cabs_pairs = generate_CABS();
  cabs_obs_ = cabs_pairs.first;
  cabs_aux_ = cabs_pairs.second;
  // Setting actual number of CABS.
  ncabs_ = cabs_obs_->mdim();

  ////////////////////
  // CABS integrals
  //////////////////

  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs(new PCompCABSFile<ERIBatch>(geom_, gamma, false, false, true, false, false, "ERI CABS"));
  eri_cabs_ = eri_cabs;
  shared_ptr<PCompCABSFile<SlaterBatch> >
    stg_cabs(new PCompCABSFile<SlaterBatch>(geom_, gamma, false, false, true, false, false, "Slater CABS"));
  stg_cabs_ = stg_cabs;
  if (!use_hy2_) {
    shared_ptr<PCompCABSFile<SlaterBatch> >
      stg_cabs2(new PCompCABSFile<SlaterBatch>(geom_, gamma, false, false, true, true, false, "Slater CABS"));
    stg_cabs2_ = stg_cabs2;
  }

#ifndef ONLY_V_X
  fill_in_cabs_matices();
#endif

  // AO ERI has been computed in the SCF class.

  // Fully transform aa/ii integrals and dump them to disk (... focus is on MP2-R12).
#ifndef ONLY_B
#ifndef ONLY_P
#ifndef ONLY_X
  eri_ii_pp_ = eri_obs_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                      nfrc_, nocc_, nfrc_, nocc_,
                                      0, nbasis_, 0, nbasis_, "ERI (pp/ii)");

  ///////////////////////////////////////////////
  // Compute the conventional MP2 contribution
  /////////////////////////////////////////////
  compute_conv_mp2();
#endif
#endif
#endif

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
#ifndef ONLY_B
#ifndef ONLY_P
    yp_ = yp;
#endif
#endif
  }


  //////////////////////////////////////
  // Yukawa ii/ii for V intermediate
  ////////////////////////////////////
#ifndef ONLY_T_Q
#ifndef ONLY_P2
#ifndef ONLY_P3
  {
#ifndef ONLY_X
#ifndef ONLY_P
    RefMOFile yp_ii_ii = yp_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                           nfrc_, nocc_, nfrc_, nocc_,
                                           nfrc_, nocc_, nfrc_, nocc_, "Yukawa (ii/ii)");
    yp_ii_ii_ = yp_ii_ii;
#endif
#endif
#ifndef SKIP_CASE1
    RefMOFile stg_ii_pp = stg_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              0, nbasis_, 0, nbasis_, "Slater (pp/ii)");
    stg_ii_pp_ = stg_ii_pp;
#endif
  }
#endif
#endif
#endif

#ifndef ONLY_B

  ///////////////////////////////
  // MO transformation for ERI
  /////////////////////////////
#ifndef ONLY_P
#ifndef ONLY_X
  {
    RefMOFile eri_ii_Ai = eri_obs_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                                  nfrc_, nocc_, nfrc_, nocc_,
                                                  0, ncabs_, 0, nocc_, "v^ia'_ii, OBS part");

    *eri_ii_Ai += *eri_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, nocc_, "v^ia'_ii, auxiliary functions");
    eri_ii_Ai_ = eri_ii_Ai;
  }
#endif
#endif

  ///////////////////////////////
  // MO transformation for STG
  /////////////////////////////
#ifndef ONLY_P
  {
    RefMOFile stg_ii_Ai = stg_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                             nfrc_, nocc_, nfrc_, nocc_,
                                             0, ncabs_, 0, nocc_, "F^ia'_ii, OBS part");
    *stg_ii_Ai += *stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, nocc_, "F^ia'_ii, auxiliary functions");
    stg_ii_Ai_ = stg_ii_Ai;
  }
#endif


  //////////////////////
  // V intermediate
  ////////////////////
#ifndef ONLY_P
#ifndef ONLY_X
  {
    RefMOFile vF = stg_ii_pp_->contract(eri_ii_pp_, "F * v (ii/ii) OBS");
    RefMOFile V_obs(new PMOFile<complex<double> >(*yp_ii_ii_ - *vF));

    RefMOFile V_cabs = stg_ii_Ai_->contract(eri_ii_Ai_, "F * v (ii/ii) CABS");

    V_cabs->flip_symmetry();
    RefMOFile V_pre(new PMOFile<complex<double> >(*V_obs - *V_cabs));

    V_ = V_pre;
  }
#endif
#endif
#endif // only_B


  /////////////////////////////////////
  // Slater & Yukawa ints with 2gamma
  ///////////////////////////////////
#ifndef ONLY_P
  {
    shared_ptr<PairCompFile<SlaterBatch> >
      stg_yp2(new PairCompFile<SlaterBatch>(geom_, 2.0 * gamma, "Slater and Yukawa ints (2gamma)"));
    stg_yp2->store_integrals();
    stg_yp2->reopen_with_inout();
    shared_ptr<PCompFile<SlaterBatch> > stg2 = stg_yp2->first();
    shared_ptr<PCompFile<SlaterBatch> > yp2  = stg_yp2->second();  // TODO delete YP from here.
    stg2_ = stg2;
  }
#endif

  ////////////////////
  // X intermediate
  //////////////////
#ifndef ONLY_P
  {
    RefMOFile stg2_ii_ii = stg2_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                               nfrc_, nocc_, nfrc_, nocc_,
                                               nfrc_, nocc_, nfrc_, nocc_, "Slater (ii/ii) 2gamma");
    stg2_ii_ii_ = stg2_ii_ii;
#ifndef ONLY_B
#ifdef ONLY_X
    {
      const complex<double> en_xtt = stg2_ii_ii_->get_energy_two_amp_X(eig_);
      cout << "**** debug ****  X(F2) contrib. " << setprecision(10) << en_xtt << endl;
    }
#endif
    RefMOFile FF = stg_ii_pp_->contract(stg_ii_pp_, "F * F (ii/ii) OBS");
    RefMOFile X_obs(new PMOFile<complex<double> >(*stg2_ii_ii - *FF));
#ifdef ONLY_X
    {
      const complex<double> en_xtt = X_obs->get_energy_two_amp_X(eig_);
      cout << "**** debug ****  X(obs) contrib. " << setprecision(10) << en_xtt << endl;
    }
#endif

    RefMOFile X_cabs = stg_ii_Ai_->contract(stg_ii_Ai_, "F * F (ii/ii) CABS");

    X_cabs->flip_symmetry();
#ifdef ONLY_X
    {
      const complex<double> en_xtt = X_cabs->get_energy_two_amp_X(eig_);
      cout << "**** debug ****  X(cabs) contrib. " << setprecision(10) << en_xtt << endl;
    }
#endif
    RefMOFile X_pre(new PMOFile<complex<double> >(*X_obs - *X_cabs));

    X_ = X_pre;
    cout << "**** debug ****  X contrib. for debug" << setprecision(10) << X_->get_energy_two_amp_B().real() << endl;
#endif // ifndef ONLY_B
  }

#ifndef ONLY_B
#ifdef DEBUG_PRINT
  {
#ifndef ONLY_X
    const complex<double> en_vt = V_->get_energy_one_amp();
    cout << "**** debug ****  V contrib. " << setprecision(10) << en_vt.real() << endl;
#endif
    const complex<double> en_xtt = X_->get_energy_two_amp_X(eig_);
    cout << "**** debug ****  X contrib. " << setprecision(10) << en_xtt.real() << endl;
  }
#endif

#endif // ifndef ONLY_B
#endif // ifndef ONLY_P

#ifndef ONLY_V_X
  /////////////////////
  // B intermediate
  ///////////////////
  // Approximation C: Kedzuch et al. IJQC 105, 929 (2005).
  {
#ifndef ONLY_P
    // T intermediate (direct)
    RefMOFile T(new PMOFile<complex<double> >(*stg2_ii_ii_ * (gamma*gamma)));

#ifdef DEBUG_PRINT
    cout << "**** debug ****  T contrib. " << setprecision(10) << T->get_energy_two_amp_B().real() << endl;
#endif

    // Q intermediate (made of X * h)
    RefMOFile Q;
    {
      RefMatrix hj_ip(new PMatrix1e(hJ_obs_obs_, make_pair(0, nocc_)));
      RefPCoeff chj_ip(new PCoeff(*coeff_ * *hj_ip));

      RefMatrix hj_iA(new PMatrix1e(hJ_obs_cabs_, make_pair(0, nocc_)));
      RefPCoeff chj_iA_comb(new PCoeff(*coeff_cabs_ * *hj_iA));

      pair<RefMatrix, RefMatrix> chj_iA = chj_iA_comb->split(geom_->nbasis(), geom_->naux());
      RefMatrix chj_iA_obs = chj_iA.first;
      RefCoeff chj_iA_cabs(new PCoeff(*chj_iA.second));

//#define ORB_DEBUG
#ifndef ORB_DEBUG
      *chj_ip += *chj_iA_obs;
#endif
      chj_iA_cabs->scale(0.5);
      chj_ip->scale(0.5);

      // MO transform using Hartree-weighted index
      RefMOFile X_ii_hi = stg2_->mo_transform(coeff_, coeff_, chj_ip, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              nfrc_, nocc_, nfrc_, nocc_, "Q intermediate: stg2 (OBS) 1/2");

      // integral evaluation...
      RefMOFile X_ii_hi_cabs;
      {
        shared_ptr<PCompCABSFile<SlaterBatch> >
          stg2_cabs(new PCompCABSFile<SlaterBatch>(geom_, 2.0 * gamma, false, false, true, false,
                                                   false, "Slater CABS (2gamma)"));
        // Use integral-direct mo transform.
        X_ii_hi_cabs = stg2_cabs->mo_transform_cabs_aux(coeff_, coeff_, chj_iA_cabs, coeff_,
                                                        nfrc_, nocc_, nfrc_, nocc_,
                                                        nfrc_, nocc_, nfrc_, nocc_, "Q intermediate: stg2, direct (CABS) 2/2",
                                                        true);
      }
#ifndef ORB_DEBUG
      *X_ii_hi += *X_ii_hi_cabs;
#endif

      X_ii_hi->flip_symmetry();
      RefMOFile Qtmp(new PMOFile<complex<double> >(*X_ii_hi));
      Q = Qtmp;
      Q->scale(2.0);
    } // end of Q intermediate construction.

#ifdef DEBUG_PRINT
    cout << "**** debug ****  Q contrib.: " << setprecision(10) << Q->get_energy_two_amp_B().real() << endl;
#endif

#endif // only-p

#ifndef ONLY_T_Q

    // P1 intermediate R^PQ_ij K^R_P R^kl_RQ
    RefMOFile P;

    {
      RefCoeff cfobs(new PCoeff(*coeff_ * *K_obs_obs_));
      RefCoeff cfri(new PCoeff(*coeff_cabs_ * *K_obs_cabs_));
      pair<RefCoeff, RefCoeff> p = cfri->split(geom_->nbasis(), geom_->naux());
      RefCoeff cfcabs_obs = p.first;
      RefCoeff cfcabs_aux = p.second;
      RefCoeff cfcabs_obs_fold(new PCoeff(*cfobs + *cfcabs_obs));

      // OBS * CABS quantity
      RefCoeff xfobs(new PCoeff(*coeff_ * *K_cabs_obs_));

      // RI * CABS quantity.
      RefCoeff xfri(new PCoeff(*coeff_cabs_ * *K_cabs_cabs_));
      // Split into OBS*CABS & CABS*CABS quantities.
      pair<RefCoeff, RefCoeff> q = xfri->split(geom_->nbasis(), geom_->naux());
      RefCoeff xfcabs_obs = q.first;
      RefCoeff xfcabs_aux = q.second;
      RefCoeff xfcabs_obs_fold(new PCoeff(*xfobs + *xfcabs_obs));

#ifndef ONLY_P2
#ifndef ONLY_P3
      RefMOFile p1;
      // first, evaluate R^Pq_ij K^r_P R^kl_rq
      {
#if 0
        RefMOFile p1_1 = stg_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, nbasis_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case1, OBS 1/3");
#endif
#ifndef SKIP_CASE1
        RefMOFile p1_2 = stg_->mo_transform(coeff_, coeff_, cfcabs_obs_fold, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, nbasis_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case1, OBS 1/2");
        *p1_2 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cfcabs_aux, coeff_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, nbasis_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case1, CABS 2/2"));
        p1 = stg_ii_pp_->contract(p1_2, "P1: R^PQ_ij K^R_P R^kl_RQ case1");
#endif
      }
#ifdef DEBUG_PRINT
#ifndef SKIP_CASE1
      cout << "**** debug ****  P1 case1 contrib.: " << setprecision(10) << p1->get_energy_two_amp_B().real() << endl;
#endif
#endif
      // second, evaluate R^PA_ij K^r_P R^kl_rA
      RefMOFile p1_5_ket;
      if (!use_hy2_) {
        RefMOFile p1_3 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, OBS 1/6");
        *p1_3 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, CABS 2/6"));
        p1_5_ket = p1_3;
        RefMOFile p1_4 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, cfcabs_obs_fold,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, OBS 3/6");
        *p1_4 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cfcabs_obs_fold,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, CABS 4/6"));
        *p1_4 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cfcabs_aux, cabs_obs_,
                                                     nfrc_, nocc_, nfrc_, nocc_,
                                                     0, nbasis_, 0, ncabs_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, CABS 5/6"
                                                     )->flip_sort());
        *p1_4 += *(stg_cabs2_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cfcabs_aux,
                                                     nfrc_, nocc_, nfrc_, nocc_,
                                                     0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, CABS2 6/6"));
#ifndef SKIP_CASE1
        *p1 += *(p1_3->contract(p1_4, "P1: R^PQ_ij K^R_P R^kl_RQ case2"));
#else
        p1 = p1_3->contract(p1_4, "P1: R^PQ_ij K^R_P R^kl_RQ case2");
#endif

      } else {
        // In HY2, R^PA_ij K^r_P R^kl_rA has been approximated by R^pA_ij K^r_p R^kl_rA
        RefMOFile p1_3 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, OBS 1/4");
        *p1_3 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, CABS 2/4"));
        p1_5_ket = p1_3;
        RefMOFile p1_4 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, cfobs,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, OBS 3/4");
        *p1_4 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cfobs,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, CABS 4/4"));
#ifndef SKIP_CASE1
        *p1 += *(p1_3->contract(p1_4, "P1: R^PQ_ij K^R_P R^kl_RQ case2"));
#else
        p1 = p1_3->contract(p1_4, "P1: R^PQ_ij K^R_P R^kl_RQ case2");
#endif
      }
#ifdef DEBUG_PRINT
      cout << "**** debug ****  P1 case2 summed: " << setprecision(10) << p1->get_energy_two_amp_B().real() << endl;
#endif
      // third, evaluate R^Pq_ij K^A_P R^kl_Aq
      {
        RefMOFile p1_5 = p1_5_ket;
        RefMOFile p1_6 = stg_->mo_transform(coeff_, coeff_, xfcabs_obs_fold, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case3, OBS 1/2");
        *p1_6 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, xfcabs_aux, coeff_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case3, CABS 2/2"));
        *p1 += *(p1_5->contract(p1_6, "P1: R^PQ_ij K^R_P R^kl_RQ case3"));
      }
#ifdef DEBUG_PRINT
      cout << "**** debug ****  P1 case3 summed: " << setprecision(10) << p1->get_energy_two_amp_B().real() << endl;
#endif
      // then, evaluate R^PB_ij K^A_P R^kl_AB
      if (!use_hy2_) {

        RefMOFile p1_7 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, cabs_obs_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, ncabs_, "P1: R^PQ_ij K^R_P R^kl_RQ case4, OBS 1/7");
        RefMOFile p1_7p = stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cabs_obs_,
                                                           nfrc_, nocc_, nfrc_, nocc_,
                                                           0, ncabs_, 0, ncabs_, "P1: R^PQ_ij K^R_P R^kl_RQ case4, CABS 2/7");
        *p1_7 += *p1_7p;
        *p1_7 += *p1_7p->flip_sort();
        *p1_7 += *(stg_cabs2_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cabs_aux_,
                                                     nfrc_, nocc_, nfrc_, nocc_,
                                                     0, ncabs_, 0, ncabs_, "P1: R^PQ_ij K^R_P R^kl_RQ case4, CABS2 3/7"));

        RefMOFile p1_8 = stg_->mo_transform(coeff_, coeff_, xfcabs_obs_fold, cabs_obs_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, ncabs_, "P1: R^PQ_ij K^R_P R^kl_RQ case4, OBS 4/7");
        *p1_8 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, xfcabs_aux, cabs_obs_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, ncabs_, "P1: R^PQ_ij K^R_P R^kl_RQ case4, CABS 5/7"));
        {
          RefMOFile p1_8p = stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, xfcabs_obs_fold,
                                                             nfrc_, nocc_, nfrc_, nocc_,
                                                             0, ncabs_, 0, ncabs_, "P1: R^PQ_ij K^R_P R^kl_RQ case4, CABS 6/7");
          *p1_8 += *(p1_8p->flip_sort());
        }
        *p1_8 += *(stg_cabs2_->mo_transform_cabs_aux(coeff_, coeff_, xfcabs_aux, cabs_aux_,
                                                     nfrc_, nocc_, nfrc_, nocc_,
                                                     0, ncabs_, 0, ncabs_, "P1: R^PQ_ij K^R_P R^kl_RQ case4, CABS2 7/7"));

        *p1 += *(p1_7->contract(p1_8, "P1: R^PQ_ij K^R_P R^kl_RQ case4"));
      }
#ifdef DEBUG_PRINT
    cout << "**** debug ****  P1 case4 summed: " << setprecision(10) << p1->get_energy_two_amp_B().real() << endl;
#endif

      p1->flip_symmetry();
      P = p1;
#ifdef DEBUG_PRINT
      cout << "**** debug ****  P1 contrib.: " << setprecision(10) << p1->get_energy_two_amp_B().real() << endl;
#endif
#endif // ONLY_P3
#endif // ONLY_P2
    }

#ifndef ONLY_P1
    RefCoeff cfobs(new PCoeff(*coeff_ * *fock_obs_obs_));

    RefCoeff cfri(new PCoeff(*coeff_cabs_ * *fock_obs_cabs_));
    pair<RefCoeff, RefCoeff> p = cfri->split(geom_->nbasis(), geom_->naux());
    RefCoeff cfcabs_obs = p.first;
    RefCoeff cfcabs_aux = p.second;
    RefCoeff cfcabs_obs_fold(new PCoeff(*cfobs + *cfcabs_obs));

    // OBS * CABS quantity
    RefCoeff xfobs(new PCoeff(*coeff_ * *fock_cabs_obs_));
    // RI * CABS quantity.
    RefCoeff xfri(new PCoeff(*coeff_cabs_ * *fock_cabs_cabs_));

    // Split into OBS*CABS & CABS*CABS quantities.
    pair<RefCoeff, RefCoeff> q = xfri->split(geom_->nbasis(), geom_->naux());
    RefCoeff xfcabs_obs = q.first;
    RefCoeff xfcabs_aux = q.second;
    RefCoeff xfcabs_obs_fold(new PCoeff(*xfobs + *xfcabs_obs));

    // P2 intermediate R^Pm_ij f^Q_P R^kl_Qm
    {
#ifndef ONLY_P3
      RefMOFile p2;
//#define USE_OLD
#ifdef USE_OLD
      // first, evaluate R^pm_ij f^Q_p R^kl_Qm
      {
        RefMOFile p2_1 = stg_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, nbasis_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=p, OBS 1/3");
        RefMOFile p2_2 = stg_->mo_transform(coeff_, coeff_, cfcabs_obs_fold, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, nbasis_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=p, OBS 2/3");
        *p2_2 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cfcabs_aux, coeff_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, nbasis_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=p, CABS 2/3"));
        p2 = p2_1->contract(p2_2, "P2: R^Pm_ij f^Q_P R^kl_Qm (P=p)");
      }
      // then, evaluate R^Am_ij f^Q_A R^kl_Qm
      {
        RefMOFile p2_1 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=A, OBS 1/4");
        *p2_1 += *stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                   nfrc_, nocc_, nfrc_, nocc_,
                                                   0, ncabs_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=A, CABS 2/4");
        RefMOFile p2_2 = stg_->mo_transform(coeff_, coeff_, xfcabs_obs_fold, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=A, OBS 3/4");
        *p2_2 += *stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, xfcabs_aux, coeff_,
                                                   nfrc_, nocc_, nfrc_, nocc_,
                                                   0, ncabs_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=A, CABS 4/4");
        *p2 += *(p2_1->contract(p2_2, "P2: R^Pm_ij f^Q_P R^kl_Qm (P=p)"));
      }
#else
      // first, evaluate R^p(P)m_ij f^Q_p(P) R^kl_Qm
      {
        const tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> ri4 = generate_RI();
        // nbasis * nbasis size
        RefCoeff ri_obs_obs(new PCoeff(*get<0>(ri4)));
        // nbasis * naux size
        RefCoeff ri_obs_cabs(new PCoeff(*get<1>(ri4)));
        // naux * nbasis size
        RefCoeff ri_cabs_obs(new PCoeff(*get<2>(ri4)));
        // naux * naux size
        RefCoeff ri_cabs_cabs(new PCoeff(*get<3>(ri4)));

        const tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> fri4 = generate_fock_weighted_RI();
        // nbasis * nbasis size
        RefCoeff fri_obs_obs(new PCoeff(*get<0>(fri4)));
        // nbasis * naux size
        RefCoeff fri_obs_cabs(new PCoeff(*get<1>(fri4)));
        // naux * nbasis size
        RefCoeff fri_cabs_obs(new PCoeff(*get<2>(fri4)));
        // naux * naux size
        RefCoeff fri_cabs_cabs(new PCoeff(*get<3>(fri4)));
        {
          // first target is OBS space
          RefMOFile p2_1 = stg_->mo_transform(coeff_, coeff_, ri_obs_obs, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              0, nbasis_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=p, OBS 1/8");
          *p2_1 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, ri_obs_cabs, coeff_,
                                                      nfrc_, nocc_, nfrc_, nocc_,
                                                      0, nbasis_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=p, CABS 2/8"));
          RefMOFile p2_2 = stg_->mo_transform(coeff_, coeff_, fri_obs_obs, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              0, nbasis_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=p, OBS 3/8");
          *p2_2 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, fri_obs_cabs, coeff_,
                                                      nfrc_, nocc_, nfrc_, nocc_,
                                                      0, nbasis_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=p, CABS 4/8"));
          p2 = p2_2->contract(p2_1, "P2: R^Pm_ij f^Q_P R^kl_Qm (P=p)");
        }
        // then, evaluate R^A(P)m_ij f^Q_A(P) R^kl_Qm
        {
          RefMOFile p2_1 = stg_->mo_transform(coeff_, coeff_, ri_cabs_obs, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              0, ncabs_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=A, OBS 5/8");
          *p2_1 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, ri_cabs_cabs, coeff_,
                                                      nfrc_, nocc_, nfrc_, nocc_,
                                                      0, ncabs_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=A, CABS 6/8"));
          RefMOFile p2_2 = stg_->mo_transform(coeff_, coeff_, fri_cabs_obs, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              0, ncabs_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=A, OBS 7/8");
          *p2_2 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, fri_cabs_cabs, coeff_,
                                                      nfrc_, nocc_, nfrc_, nocc_,
                                                      0, ncabs_, 0, nocc_, "P2: R^Pm_ij f^Q_P R^kl_Qm, P=A, CABS 8/8"));
          *p2 += *(p2_2->contract(p2_1, "P2: R^Pm_ij f^Q_P R^kl_Qm (P=p)"));
        }
      }
#endif
      p2->flip_symmetry();
#ifndef ONLY_P2
      *P += *p2;
#endif

#ifdef DEBUG_PRINT
      cout << "**** debug ****  P2 contrib.: " << setprecision(10) << p2->get_energy_two_amp_B().real() << endl;
#endif
#endif // ONLY_P3
    }

#ifndef ONLY_P2
    // P3 intermediate R^mA_ij f^n_m R^kl_nA;
    {
      RefMOFile p3_ket;
      {
        RefMOFile p3_1 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, cfobs,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, nocc_, "P3: R^Am_ij f^m_n R^kl_An, 1/2, OBS");
        *p3_1 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cfobs,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, nocc_, "P3: R^mA_ij f^m_n R^kl_nA, 1/2, CABS"));
        RefMOFile p3_2 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, nocc_, "P3: R^Am_ij f^m_n R^kl_An, 2/2, OBS");
        *p3_2 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, nocc_, "P3: R^mA_ij f^m_n R^kl_nA, 2/2, CABS"));
        p3_ket = p3_2;
        RefMOFile p3 = p3_1->contract(p3_2, "P3: R^Am_ij f^m_n R^kl_An");

        p3->flip_symmetry();

#ifndef ONLY_P3
        if (use_hy2_) {
          // in HY2, -P3 + 2*P5A is approximated by +P3, which is quite accurate.
          *P += *p3;
        } else {
          *P -= *p3;
        }
#endif
#ifdef DEBUG_PRINT
        cout << "**** debug ****  P3 contrib.: " << setprecision(10) << p3->get_energy_two_amp_B().real() << endl;
#endif
      }

       // P5A intermediate R^mA_ij f^P_m R^kl_PA
      if (!use_hy2_) {
        RefMOFile p5a_1 = stg_cabs2_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cfcabs_aux,
                                                            nfrc_, nocc_, nfrc_, nocc_,
                                                            0, ncabs_, 0, nocc_, "P5A: R^mA_ij f^P_m R^kl_PA, 1/4 CABS_CABS");
        *p5a_1 += *(stg_->mo_transform(coeff_, coeff_, cabs_obs_, cfcabs_obs_fold,
                                       nfrc_, nocc_, nfrc_, nocc_,
                                       0, ncabs_, 0, nocc_, "P5A: R^mA_ij f^P_m R^kl_PA, 2/4 OBS_OBS"));
        *p5a_1 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, cfcabs_obs_fold,
                                                     nfrc_, nocc_, nfrc_, nocc_,
                                                     0, ncabs_, 0, nocc_, "P5A: R^mA_ij f^P_m R^kl_PA, 3/4 CABS_OBS"));
        RefMOFile tmp = stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cfcabs_aux, cabs_obs_,
                                                         nfrc_, nocc_, nfrc_, nocc_,
                                                         0, nocc_, 0, ncabs_, "P5A: R^mA_ij f^P_m R^kl_PA, 4/4 OBS_CABS");
        *p5a_1 += *(tmp->flip_sort());
        RefMOFile p5a = p5a_1->contract(p3_ket, "P5A: R^mA_ij f^P_m R^kl_PA");
        p5a->flip_symmetry();
        p5a->scale(2.0);
        *P += *p5a;
#ifdef DEBUG_PRINT
        cout << "**** debug ****  P5a contrib.: " << setprecision(10) << p5a->get_energy_two_amp_B().real() << endl;
#endif
      }
    }


#ifndef ONLY_P3
    // P4 intermediate R^pb_ij f^r_p R^kl_rb
    {
      RefMOFile p4_ket;
      {
        RefMOFile p4_1 = stg_->mo_transform(coeff_, coeff_, cfobs, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, nbasis_, nocc_ , nbasis_, "P4: R^pb_ij f^r_p R^kl_rb, 1/2");
        RefMOFile p4_2 = stg_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, nbasis_, nocc_ , nbasis_, "P4: R^pb_ij f^r_p R^kl_rb, 2/2");
        if (true) p4_ket = p4_2;
        //if (!use_hy2_) p4_ket = p4_2;

        RefMOFile p4 = p4_1->contract(p4_2, "P4: R^pb_ij f^r_p R^ij_rb");
        p4->flip_symmetry();
        *P += *p4;
#ifdef DEBUG_PRINT
        cout << "**** debug ****  P4 contrib.: " << setprecision(10) << p4->get_energy_two_amp_B().real() << endl;
#endif
      }

      // P5b intermediate  R^Ab_ij f^p_A R^kl_pb
      {
        RefMOFile p5b;
        //if (!use_hy2_) {
        if (true) {

          RefMOFile p5b_1 = stg_->mo_transform(coeff_, coeff_, cfcabs_obs, coeff_,
                                               nfrc_, nocc_, nfrc_, nocc_,
                                               0, nbasis_, nocc_, nbasis_, "P5B: R^Ab_ij f^p_A R^kl_pb, 1/2, OBS");
          *p5b_1 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cfcabs_aux, coeff_,
                                                       nfrc_, nocc_, nfrc_, nocc_,
                                                       0, nbasis_, nocc_, nbasis_, "P5B: R^Ab_ij f^p_A R^kl_pb, 2/2, CABS"));
          p5b = p5b_1->contract(p4_ket, "R^Ab_ij f^p_A R^kl_pb");

        } else {
          // P5b intermediate  R^Ab_ij f^p_A R^kl_pb approximated by R^Ab_ij f^a_A R^kl_ab
          RefMOFile p5b_1 = stg_->mo_transform(coeff_, coeff_, cfcabs_obs, coeff_,
                                               nfrc_, nocc_, nfrc_, nocc_,
                                               nocc_, nbasis_, nocc_, nbasis_, "P5B: R^Ab_ij f^p_A R^kl_pb, 1/3, OBS");
          *p5b_1 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cfcabs_aux, coeff_,
                                                       nfrc_, nocc_, nfrc_, nocc_,
                                                       nocc_, nbasis_, nocc_, nbasis_, "P5B: R^Ab_ij f^p_A R^kl_pb, 2/3, CABS"));
          RefMOFile p5b_2 = stg_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                               nfrc_, nocc_, nfrc_, nocc_,
                                               nocc_, nbasis_, nocc_, nbasis_, "P5B: R^Ab_ij f^p_A R^kl_pb, 3/3, OBS");
          p5b = p5b_1->contract(p5b_2, "R^Ab_ij f^p_A R^kl_pb");
        }
        p5b->flip_symmetry();
        p5b->scale(2.0);
        *P += *p5b;

#ifdef DEBUG_PRINT
        cout << "**** debug ****  P5b contrib.: " << setprecision(10) << p5b->get_energy_two_amp_B().real() << endl;
#endif
      }
    }

#ifndef ONLY_P
    RefMOFile btmp(new PMOFile<complex<double> >(*T + *Q - *P));
    B_ = btmp;
#endif

#endif // ONLY_P3
#endif // ONLY_P2
#endif
#endif // ifndef ONLY_T_Q

  } // end of B intermediate construction.


  ///////////////////////////////////////
  // R12 energy contribution
  /////////////////////////////////////
#ifndef ONLY_B
#ifndef ONLY_P
  const complex<double> en_vt = V_->get_energy_one_amp();
  const complex<double> en_xtt = X_->get_energy_two_amp_X(eig_);
  const complex<double> en_btt = B_->get_energy_two_amp_B();
  cout << "  ------------------------------------------" << endl;
  cout << "    R12 contribution  : " << fixed << setw(15) << setprecision(10) << (en_vt * 2.0 + en_btt - en_xtt).real()
       << endl << endl;
  cout << "      >> V*t          : " << fixed << setw(15) << setprecision(10) << en_vt.real() << endl;
  cout << "      >> t*(ei+ej)X*t : " << fixed << setw(15) << setprecision(10) << en_xtt.real() << endl;
  cout << "      >> t*B*t        : " << fixed << setw(15) << setprecision(10) << en_btt.real() << endl;
  cout << "  ------------------------------------------" << endl << endl;

//#define LOCAL_DEBUG_PMP2
#ifdef LOCAL_DEBUG_PMP2
  V_->rprint();
  X_->rprint();
  B_->rprint();
#endif
#endif // ifndef ONLY_P2
#endif
#endif // ifndef ONLY_B
}


