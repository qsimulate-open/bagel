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
    RefMOFile yp_ii_ii = yp_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            nfrc_, nocc_, nfrc_, nocc_, "Yukawa (ii/ii)");
    yp_ii_ii_ = yp_ii_ii;
    RefMOFile stg_ii_pp = stg_->mo_transform(coeff_, coeff_, coeff_, coeff_,
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
  {
    shared_ptr<PCompCABSFile<SlaterBatch> >
      stg_cabs2(new PCompCABSFile<SlaterBatch>(geom_, gamma, false, false, true, true, false, "Slater CABS"));
    stg_cabs2->store_integrals();
    stg_cabs2->reopen_with_inout();
    stg_cabs2_ = stg_cabs2;
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
    RefMOFile eri_ii_pi = eri_obs_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                                  nfrc_, nocc_, nfrc_, nocc_,
                                                  0, ncabs_, 0, nocc_, "v^ia'_ii, OBS part");

    RefMOFile eri_ii_xi = eri_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                            nfrc_, nocc_, nfrc_, nocc_,
                                                            0, ncabs_, 0, nocc_, "v^ia'_ii, auxiliary functions");
    RefMOFile eri_ii_Ai(new PMOFile<complex<double> >(*eri_ii_xi + *eri_ii_pi));
    eri_ii_Ai_ = eri_ii_Ai;
  }

  ///////////////////////////////
  // MO transformation for STG
  /////////////////////////////
  {
    RefMOFile stg_ii_pi = stg_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              0, ncabs_, 0, nocc_, "F^ia'_ii, OBS part");
    RefMOFile stg_ii_xi = stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                            nfrc_, nocc_, nfrc_, nocc_,
                                                            0, ncabs_, 0, nocc_, "F^ia'_ii, auxiliary functions");
    RefMOFile stg_ii_Ai(new PMOFile<complex<double> >(*stg_ii_xi + *stg_ii_pi));
    stg_ii_Ai_ = stg_ii_Ai;
  }


  //////////////////////
  // V intermediate
  ////////////////////
  {
    RefMOFile vF = stg_ii_pp_->contract(eri_ii_pp_, "F * v (ii/ii) OBS");
    RefMOFile V_obs(new PMOFile<complex<double> >(*yp_ii_ii_ - *vF));

    RefMOFile V_cabs = stg_ii_Ai_->contract(eri_ii_Ai_, "F * v (ii/ii) CABS");

    V_cabs->flip_symmetry();
    RefMOFile V_pre(new PMOFile<complex<double> >(*V_obs - *V_cabs));

    V_ = V_pre;
  }


  ///////////////////////////////////////
  // Direct contribution to R12 energy
  /////////////////////////////////////
  cout << "  * Built V intermediate." << endl;
  complex<double> en_vt = V_->get_energy_one_amp();
  cout << "  F12 energy (Vt): " << setprecision(10) << en_vt.real() << endl << endl;


  /////////////////////////////////////
  // Slater & Yukawa ints with 2gamma
  ///////////////////////////////////
  {
    shared_ptr<PairCompFile<SlaterBatch> >
      stg_yp2(new PairCompFile<SlaterBatch>(geom_, 2.0 * gamma, "Slater and Yukawa ints (2gamma)"));
    stg_yp2->store_integrals();
    stg_yp2->reopen_with_inout();
    shared_ptr<PCompFile<SlaterBatch> > stg2 = stg_yp2->first();
    shared_ptr<PCompFile<SlaterBatch> > yp2  = stg_yp2->second();
    stg2_ = stg2;
  }

  ////////////////////
  // X intermediate
  //////////////////
  {
    RefMOFile stg2_ii_ii = stg2_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                               nfrc_, nocc_, nfrc_, nocc_,
                                               nfrc_, nocc_, nfrc_, nocc_, "Slater (ii/ii) 2gamma");
    stg2_ii_ii_ = stg2_ii_ii;
    RefMOFile FF = stg_ii_pp_->contract(stg_ii_pp_, "F * F (ii/ii) OBS");
    RefMOFile X_obs(new PMOFile<complex<double> >(*stg2_ii_ii - *FF));

    RefMOFile X_cabs = stg_ii_Ai_->contract(stg_ii_Ai_, "F * F (ii/ii) CABS");

    X_cabs->flip_symmetry();
    RefMOFile X_pre(new PMOFile<complex<double> >(*X_obs - *X_cabs));

    X_ = X_pre;
  }


  /////////////////////
  // B intermediate
  ///////////////////
  // Approximation C: Kedzuch et al. IJQC 105, 929 (2005).
  {
    // T intermediate (direct)
    RefMOFile T(new PMOFile<complex<double> >(*stg2_ii_ii_ * (gamma*gamma)));

    // Q intermediate (made of X * h)
    RefMOFile Q;
    {
      // Hartree builder (needs modification!!!)
      // nbasis * nbasis size
      hJ_obs_obs_ = generate_hJ_obs_obs();
      // nbasis * ncabs size
      hJ_obs_cabs_ = generate_hJ_obs_cabs();

      pair<RefMatrix, RefMatrix> p = generate_hJ_cabs_pair();
      // ncabs * nbasis size
      hJ_cabs_obs_ = p.first;
      // ncabs * ncabs size
      hJ_cabs_cabs_ = p.second;


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
      RefMOFile X_ii_ih = stg2_->mo_transform(coeff_, coeff_,  chj_ip, coeff_,
                                              nfrc_, nocc_, nfrc_, nocc_,
                                              nfrc_, nocc_, nfrc_, nocc_, "Q intermediate: stg2 (OBS) 1/2");

      // integral evaluation...
      RefMOFile X_ii_ih_cabs;
      {
        shared_ptr<PCompCABSFile<SlaterBatch> >
          stg2_cabs(new PCompCABSFile<SlaterBatch>(geom_, 2.0 * gamma, false, false, true, false,
                                                   false, "Slater CABS (2gamma)"));
        stg2_cabs->store_integrals();
        stg2_cabs->reopen_with_inout();
        X_ii_ih_cabs = stg2_cabs->mo_transform_cabs_aux(coeff_, coeff_,  chj_iA_cabs, coeff_,
                                                        nfrc_, nocc_, nfrc_, nocc_,
                                                        nfrc_, nocc_, nfrc_, nocc_, "Q intermediate: stg2 (CABS) 1/2");
      }

      *X_ii_ih += *X_ii_ih_cabs;

      X_ii_ih->flip_symmetry();
      RefMOFile Qtmp(new PMOFile<complex<double> >(*X_ii_ih));
      Q = Qtmp;
      Q->scale(2.0);
    } // end of Q intermediate construction.

    //Q->rprint();

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


    // P1 intermediate R^PQ_ij K^R_P R^kl_RQ
    RefMOFile P;
    {
      RefCoeff cfobs(new PCoeff(*coeff_ * *K_obs_obs_));
      RefCoeff cfri(new PCoeff(*coeff_cabs_ * *K_obs_cabs_));
      pair<RefCoeff, RefCoeff> p = cfri->split(geom_->nbasis(), geom_->ncabs());
      RefCoeff cfcabs_obs = p.first;
      RefCoeff cfcabs_aux = p.second;
      RefCoeff cfcabs_obs_fold(new PCoeff(*cfobs + *cfcabs_obs));

      // OBS * CABS quantity
      RefCoeff xfobs(new PCoeff(*coeff_ * *K_cabs_obs_));

      // RI * CABS quantity.
      RefCoeff xfri(new PCoeff(*coeff_cabs_ * *K_cabs_cabs_));
      // Split into OBS*CABS & CABS*CABS quantities.
      pair<RefCoeff, RefCoeff> q = xfri->split(geom_->nbasis(), geom_->ncabs());
      RefCoeff xfcabs_obs = q.first;
      RefCoeff xfcabs_aux = q.second;
      RefCoeff xfcabs_obs_fold(new PCoeff(*xfobs + *xfcabs_obs));

      RefMOFile p1;
      // first, evaluate R^Pq_ij K^r_P R^kl_rq
      {
        RefMOFile p1_1 = stg_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, nbasis_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case1, OBS 1/3");
        RefMOFile p1_2 = stg_->mo_transform(coeff_, coeff_, cfcabs_obs_fold, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, nbasis_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case1, OBS 2/3");
        *p1_2 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cfcabs_aux, coeff_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, nbasis_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case1, CABS 3/3"));
        p1 = p1_1->contract(p1_2, "P1: R^PQ_ij K^R_P R^kl_RQ case1");
      }
      // second, evaluate R^PA_ij K^r_P R^kl_rA
      RefMOFile p1_5_ket;
      {
        RefMOFile p1_3 = stg_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                            nfrc_, nocc_, nfrc_, nocc_,
                                            0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, OBS 1/6");
        *p1_3 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                    nfrc_, nocc_, nfrc_, nocc_,
                                                    0, ncabs_, 0, nbasis_, "P1: R^PQ_ij K^R_P R^kl_RQ case2, CABS 2/6"));
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
        *p1 += *(p1_3->contract(p1_4, "P1: R^PQ_ij K^R_P R^kl_RQ case2"));

        p1_5_ket = p1_3;
      }
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
      // then, evaluate R^PB_ij K^A_P R^kl_AB
      {
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

        *p1 += *(p1_7->contract(p1_8, "P1: R^PQ_ij K^R_P R^kl_RQ case3"));
      }
      p1->flip_symmetry();
      P = p1;
    }


    RefCoeff cfobs(new PCoeff(*coeff_ * *fock_obs_obs_));

    RefCoeff cfri(new PCoeff(*coeff_cabs_ * *fock_obs_cabs_));
    pair<RefCoeff, RefCoeff> p = cfri->split(geom_->nbasis(), geom_->ncabs());
    RefCoeff cfcabs_obs = p.first;
    RefCoeff cfcabs_aux = p.second;
    RefCoeff cfcabs_obs_fold(new PCoeff(*cfobs + *cfcabs_obs));

    // OBS * CABS quantity
    RefCoeff xfobs(new PCoeff(*coeff_ * *fock_cabs_obs_));

    // RI * CABS quantity.
    RefCoeff xfri(new PCoeff(*coeff_cabs_ * *fock_cabs_cabs_));
    // Split into OBS*CABS & CABS*CABS quantities.
    pair<RefCoeff, RefCoeff> q = xfri->split(geom_->nbasis(), geom_->ncabs());
    RefCoeff xfcabs_obs = q.first;
    RefCoeff xfcabs_aux = q.second;
    RefCoeff xfcabs_obs_fold(new PCoeff(*xfobs + *xfcabs_obs));

    // P2 intermediate R^Pm_ij f^Q_P R^kl_Qm
    {
      RefMOFile p2;
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
      p2->flip_symmetry();
      *P += *p2;
    }

    // P3 intermediate R^mA_ij f^n_m R^kl_nA
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
        p3->scale(-1.0);

        *P += *p3;
      }

       // P5A intermediate R^mA_ij f^P_m R^kl_PA
      {
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
        //p5a->rprint();
        *P += *p5a;
      }
    }

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
        p4_ket = p4_2;
        RefMOFile p4 = p4_1->contract(p4_2, "P4: R^pb_ij f^r_p R^ij_rb");
        p4->flip_symmetry();
        //p4->rprint();

        *P += *p4;
      }

      // P5b intermediate  R^Ab_ij f^p_A R^kl_pb
      {
        RefMOFile p5b_1 = stg_->mo_transform(coeff_, coeff_, cfcabs_obs, coeff_,
                                             nfrc_, nocc_, nfrc_, nocc_,
                                             0, nbasis_, nocc_, nbasis_, "P5: R^Ab_ij f^p_A R^kl_pb, 1/1, OBS");
        *p5b_1 += *(stg_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cfcabs_aux, coeff_,
                                                     nfrc_, nocc_, nfrc_, nocc_,
                                                     0, nbasis_, nocc_, nbasis_, "R^Ab_ij f^p_A R^kl_pb, 1/1, CABS"));
        RefMOFile p5b = p5b_1->contract(p4_ket, "R^Ab_ij f^p_A R^kl_pb");
        p5b->flip_symmetry();
        p5b->scale(2.0);

        *P += *p5b;
      }
    }

    RefMOFile btmp(new PMOFile<complex<double> >(*T + *Q - *P));
    B_ = btmp;

  } // end of B intermediate construction.

#define LOCAL_DEBUG_PMP2
#ifdef LOCAL_DEBUG_PMP2
  V_->rprint();
  X_->rprint();
  B_->rprint();
#endif

}


