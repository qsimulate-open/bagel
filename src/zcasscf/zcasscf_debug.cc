//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf_debug.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/zcasscf/zcasscf.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;

static const int morbital = 1;
static const int norbital = 0;

void ZCASSCF::___debug___orbital_rotation(const bool kramers) {
  // currently transforming n and mth orbitals;
  const double angle = idata_->get<double>("debugrot", 0.0);
  if (angle == 0.0) return;

  cout << ">>>>>>>>>>>> debug >>>>>>>>>>>>" << endl;
  cout << "perturbing an orbital by : " << setprecision(5) << angle << endl;
  cout << "kramers adaptation       : " << (kramers ? "on" : "off") << endl;

  shared_ptr<ZRotFile> atmp = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2, /*superci*/false);

  // currently ++ and -- blocks
  const bool perturb_active = idata_->get<bool>("perturb_active", false);
  const bool perturb_ca     = idata_->get<bool>("perturb_ca", false);
  if (perturb_active) {
    cout << "perturbing virtual active block" << endl;
    atmp->ele_va(morbital, norbital) = angle;
    if (kramers)
      atmp->ele_va(nvirt_+morbital, nact_+norbital) = angle;
  } else if (perturb_ca) {
    cout << "perturbing closed active block" << endl;
    atmp->ele_ca(morbital, norbital) = angle;
    if (kramers)
      atmp->ele_ca(nclosed_+morbital, nact_+norbital) = angle;
  } else {
    cout << "perturbing virtual closed block" << endl;
    atmp->ele_vc(morbital, norbital) = angle;
    if (kramers)
      atmp->ele_vc(nvirt_+morbital, nclosed_+norbital) = angle; // = conj(angle)
  }

  shared_ptr<ZMatrix> amattmp = atmp->unpack<ZMatrix>();

  // multiply -1 from the formula. multiply -i to make amat hermite (will be compensated)
  *amattmp *= -1.0 * complex<double>(0.0, -1.0);

  // restore the matrix from RotFile
  unique_ptr<double[]> eig(new double[amattmp->ndim()]);
  amattmp->diagonalize(eig.get());

  auto amattmp_sav = amattmp->copy();
  for (int i = 0; i != amattmp->ndim(); ++i) {
    complex<double> ex = exp(complex<double>(0.0, eig[i]));
    for_each(amattmp->element_ptr(0,i), amattmp->element_ptr(0,i+1), [&ex](complex<double>& a) { a *= ex; });
  }

  auto expatmp = make_shared<ZMatrix>(*amattmp ^ *amattmp_sav);

  if (kramers) {
    kramers_adapt(expatmp);
    assert( ( (*expatmp->get_submatrix(nocc_*2, 0, nvirt_, nclosed_)) 
             - (*expatmp->get_submatrix(nocc_*2+nvirt_, nclosed_, nvirt_, nclosed_)->get_conjg()) ).rms() < 1e-20 );
  }

  expatmp->purify_unitary(); //FIXME ; needed to ensure orthonormality in the perturbed coefficients
  shared_ptr<const ZMatrix> pco = make_shared<const ZMatrix>(*coeff_ * *expatmp);
  auto tt = make_shared<ZMatrix>(*coeff_ - *pco);
  // check orthonormality
  {
    auto overlap = make_shared<const RelOverlap>(geom_);
    auto orth = make_shared<ZMatrix>(*pco % *overlap * *pco);
    auto unit = orth->clone();
    unit->unit();
    *orth -= *unit; 
    if (orth->rms() > 1e-14) {
      stringstream ss;
      ss << "coefficients are not orthonormal. Residual: " << scientific << setprecision(4) << orth->rms();
      throw runtime_error(ss.str());
    }
  }

  cout << "norm coeff  = "        << setprecision(20) << coeff_->norm() << endl; 
  cout << "norm pcoeff = "        << setprecision(20) << pco->norm() << endl; 
  cout << "norm delta(coeff) = "  << setprecision(8) << tt->norm() << endl;
  cout << "<<<<<<<<<<<< debug <<<<<<<<<<<<" << endl << endl;

  coeff_ = pco;
}


void ZCASSCF::___debug___print_gradient(shared_ptr<const ZRotFile> grad, const bool with_kramers) const {
  const bool perturb_active = idata_->get<bool>("perturb_active", false);
  const bool perturb_ca     = idata_->get<bool>("perturb_ca", false);
  cout << ">>>>>>>>>>>> debug >>>>>>>>>>>>" << endl;

  if (perturb_active) {
    cout << "virtual-active orbital gradient" << endl;
    // currently ++ and -- blocks
    const complex<double> gradient = grad->ele_va(morbital, norbital)
                                   + (with_kramers ? grad->ele_va(nvirt_+morbital, nact_+norbital) : complex<double>(0.0));
    cout << setprecision(10) << gradient << endl;
  } else if (perturb_ca) {
    cout << "closed-active orbital gradient" << endl;
    // currently ++ and -- blocks
    const complex<double> gradient = grad->ele_ca(morbital, norbital)
                                   + (with_kramers ? grad->ele_ca(nclosed_+morbital, nact_+norbital) : complex<double>(0.0));
    cout << setprecision(10) << gradient << endl;
  } else {
    cout << "virtual-closed orbital gradient" << endl;
    // currently ++ and -- blocks
    const complex<double> gradient = grad->ele_vc(morbital, norbital)
                                   + (with_kramers ? grad->ele_vc(nvirt_+morbital, nclosed_+norbital) : complex<double>(0.0));
    cout << setprecision(10) << gradient << endl;
  }

  cout << "<<<<<<<<<<<< debug <<<<<<<<<<<<" << endl << endl;
}


void ZCASSCF::___debug___compute_hessian(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr, const bool with_kramers) const {
  const bool perturb_active = idata_->get<bool>("perturb_active", false);
  const bool perturb_ca     = idata_->get<bool>("perturb_ca", false);
  const bool verbose        = idata_->get<bool>("verbose", false);

  shared_ptr<const ZMatrix> coeffa = coeff_->slice(nocc_*2, coeff_->mdim());
  shared_ptr<const ZMatrix> coeffi = coeff_->slice(0, nclosed_*2);
  shared_ptr<const ZMatrix> coefft = coeff_->slice(nclosed_*2, nocc_*2);

  shared_ptr<ZMatrix> cfockd;
  if (perturb_active) { // virtual-active block
    shared_ptr<const ZMatrix> rdm1 = transform_rdm1();
    cfockd = make_shared<ZMatrix>(*cfock->get_submatrix(nclosed_*2, nclosed_*2, nact_*2, nact_*2) * *rdm1);
    cfockd->hermite();

    shared_ptr<ZMatrix> maatt = ___debug___diagonal_2rdm_contraction_coulomb(coeffa);
    shared_ptr<ZMatrix> matta = ___debug___diagonal_2rdm_contraction_exchange(coeffa);
    *maatt -= *matta;

    shared_ptr<ZMatrix> mapattp = ___debug___diagonal_integrals_coulomb_active_kramers(coeffa, coefft); // 0 by symmetry
    shared_ptr<ZMatrix> maptpta = ___debug___diagonal_integrals_exchange_active_kramers(coeffa, coefft);
    *mapattp -= *maptpta; // <- need - sign for ++/-- debug ; + sign for +-/-+

    for (int t = 0; t != nact_*2; ++t) {
      for (int a = 0; a != nvirt_*2; ++a) {
        const int na = a + nocc_*2;
        (*maatt)(a, t) += ((*cfock)(na,na) * (*rdm1)(t,t))  - (*cfockd)(t,t) - (*qxr)(nclosed_*2 + t, t)
                          + (*mapattp)(a,t);
      }
    }
    *maatt += *maatt->get_conjg(); // due to ++ <-> -- and +- <-> -+ symmetry

    cout << ">>>>>>>>>>>>>>> debug >>>>>>>>>>>>>>>" << endl;
    cout << "virtual-active diagonal hessian value" << endl;
    cout << setprecision(10) << (*maatt)(morbital, norbital)*2.0 << endl;
    cout << "<<<<<<<<<<<<<<< debug <<<<<<<<<<<<<<<" << endl << endl;

  } else if (perturb_ca) { // closed-active block
    shared_ptr<const ZMatrix> rdm1 = transform_rdm1();
    cfockd = make_shared<ZMatrix>(*cfock->get_submatrix(nclosed_*2, nclosed_*2, nact_*2, nact_*2) * *rdm1);
    cfockd->hermite();

    // (1) G(1,1)_(ti,ti)
    shared_ptr<ZMatrix> mitti = ___debug___closed_active_diagonal_hessian(coeffi, coefft, cfock, afock, qxr, verbose);

    // (1.333) G(1,1)_(kt ki,ti)
    shared_ptr<ZMatrix> kmitti = ___debug___closed_active_diagonal_hessian_kramers(coeffi, coefft, verbose);

    // (1.667) G(1,2)_(kt i,kt i)
    shared_ptr<ZMatrix> kmitit = ___debug___closed_active_offdiagonal_hessian_kramers(coeffi, coefft, verbose);

    *mitti += *kmitti + *kmitit;
    *mitti += *mitti->get_conjg(); // from symmetry

    cout << ">>>>>>>>>>>>>>> debug >>>>>>>>>>>>>>" << endl;
    cout << "closed-active diagonal hessian value" << endl;
    cout << setprecision(10) << (*mitti)(morbital, norbital)*2.0 << endl;
    cout << "<<<<<<<<<<<<<<< debug <<<<<<<<<<<<<<" << endl << endl;

  } else { // virtual-closed block
    shared_ptr<ZMatrix> maaii = ___debug___diagonal_integrals_coulomb(coeffa, coeffi);
    shared_ptr<ZMatrix> maiia = ___debug___diagonal_integrals_exchange(coeffa, coeffi);
    *maiia -= *maaii;

    for (int i = 0; i != nclosed_*2; ++i) {
      for (int a = 0; a != nvirt_*2; ++a) {
        const int na = a + nocc_*2;
        (*maiia)(a, i) += (*cfock)(na,na) + (*afock)(na,na) - (*cfock)(i,i) - (*afock)(i,i);
      }
    }

    if (with_kramers) {
      shared_ptr<ZMatrix> kaaii = ___debug___diagonal_integrals_coulomb_kramers(coeffa, coeffi);
      shared_ptr<ZMatrix> kaiia = ___debug___diagonal_integrals_exchange_kramers(coeffa, coeffi);
      *maiia += *kaiia; // <- + sign for ++/-- debugging ; need - sign for +-/-+ debugging
      *maiia -= *kaaii; // this appears zero for (at least) coulomb integrals due to symmetry
      *maiia += *maiia->get_conjg(); // due to ++ <-> -- and +- <-> -+ symmetry
    }

    cout << ">>>>>>>>>>>>>>> debug >>>>>>>>>>>>>>>" << endl;
    cout << "virtual-closed diagonal hessian value" << endl;
    cout << setprecision(10) << (*maiia)(morbital, norbital)*2.0 << endl;
    cout << "<<<<<<<<<<<<<<< debug <<<<<<<<<<<<<<<" << endl << endl;
  }

  if (nact_) {
    double fcienergy = ___debug___recompute_fci_energy(cfock->get_submatrix(nclosed_*2,nclosed_*2,nact_*2,nact_*2));
    cout << ">>>>>>>>>>>> debug >>>>>>>>>>>>" << endl;
    cout << "recomputed FCI energy" << endl;
    cout << setprecision(12) << fcienergy << endl;
    cout << "<<<<<<<<<<<< debug <<<<<<<<<<<<" << endl << endl;
  }
}


double ZCASSCF::___debug___recompute_fci_energy(shared_ptr<const ZMatrix> cfock) const {
  // returns FCI energy ; requires core fock matrix for active orbitals as input
  assert(cfock->ndim() == 2*nact_ && cfock->ndim() == cfock->mdim());

  shared_ptr<const ZMatrix> rdm2 = fci_->rdm2_av();
  shared_ptr<const ZMatrix> rdm1 = transform_rdm1();
  shared_ptr<ZMatrix> ijkl = ___debug___all_integrals_coulomb_active(coeff_->slice(nclosed_*2, nclosed_*2 + nact_*2));

  double twoelen = (rdm2->dot_product(ijkl)).real();
  double oneelen = (rdm1->dot_product(cfock)).real();
  double tote = fci_->core_energy()  + geom_->nuclear_repulsion() + oneelen + 0.5*twoelen;

  return tote;
}


shared_ptr<ZMatrix> ZCASSCF::___debug___closed_active_diagonal_hessian(shared_ptr<const ZMatrix> coeffi, shared_ptr<const ZMatrix> coefft, shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr, const bool verbose) const {
  /* returns Mat(i,t) = G^{(1,1)}_{ti,ti} = cfock(tt) + afock(tt) - cfock(ii) - afock(ii) - cfockd(tt) + D(tt)*cfock(ii)
                                          - Q^{*}_{tt} + [ (ii|vu) - (iu|vi) ] G(ttvu)
                                          + ([ (ii|tu) - (iu|ti) ] D(tu) + c.c. )
                                          + [ (ti|it) - (ii|tt) ]
  for the time being, we implement it in the worst possible way... to be updated to make it efficient. */
  assert(coefft->mdim() == nact_*2);
  if (verbose)
    cout << "<<<<<<<<<<<< debug : G(1,1)(ti,ti) <<<<<<<<<<<<" << endl;

  shared_ptr<const ZMatrix> rdm1 = transform_rdm1();
  shared_ptr<ZMatrix> cfockd = make_shared<ZMatrix>(*cfock->get_submatrix(nclosed_*2, nclosed_*2, nact_*2, nact_*2) * *rdm1);
  cfockd->hermite();
  if (verbose)
    rdm1->print("1rdm");

  shared_ptr<ZMatrix> miitt = ___debug___diagonal_integrals_coulomb(coeffi, coefft);
  shared_ptr<ZMatrix> mitti = ___debug___diagonal_integrals_exchange(coeffi, coefft);
  if (verbose) {
    miitt->get_submatrix(morbital,norbital,1,1)->print("(ii|tt)");
    mitti->get_submatrix(morbital,norbital,1,1)->print("(it|ti)");
  }
  *mitti -= *miitt;

  shared_ptr<ZMatrix> miitt2rdm = ___debug___diagonal_2rdm_contraction_coulomb(coeffi);
  shared_ptr<ZMatrix> mitti2rdm = ___debug___diagonal_2rdm_contraction_exchange(coeffi);
  if (verbose) {
    miitt2rdm->get_submatrix(morbital,norbital,1,1)->print("(ii|vu) * G(ttvu)");
    mitti2rdm->get_submatrix(morbital,norbital,1,1)->print("(iu|vi) * G(ttvu)");
  }
  *miitt2rdm -= *mitti2rdm;

  shared_ptr<ZMatrix> miitt1rdm = ___debug___diagonal_1rdm_contraction_coulomb(coeffi, coefft);
  shared_ptr<ZMatrix> mitti1rdm = ___debug___diagonal_1rdm_contraction_exchange(coeffi, coefft);
  if (verbose) {
    miitt1rdm->get_submatrix(morbital,norbital,1,1)->print("(ii|tu) * D(tu)");
    mitti1rdm->get_submatrix(morbital,norbital,1,1)->print("(iu|ti) * D(tu)");
  }
  *miitt1rdm -= *mitti1rdm;
  *miitt1rdm += *miitt1rdm->get_conjg();
  if (verbose) {
    cout << setprecision(12) << "cfockdtt    = " <<  cfockd->element(norbital,norbital) << endl;
    cout << setprecision(12) << "cfockii     = " <<  cfock->element(morbital,morbital) << endl;
    cout << setprecision(12) << "afockii     = " <<  afock->element(morbital,morbital) << endl;
    cout << setprecision(12) << "cfocktt     = " <<  cfock->element(nclosed_*2+norbital,nclosed_*2+norbital) << endl;
    cout << setprecision(12) << "afocktt     = " <<  afock->element(nclosed_*2+norbital,nclosed_*2+norbital) << endl;
    cout << setprecision(12) << "rdm1tt      = " <<  rdm1->element(norbital,norbital) << endl;
    cout << setprecision(12) << "qvectt      = " <<  qxr->element(nclosed_*2 + norbital, norbital) << endl;
  }

  for (int t = 0; t != nact_*2; ++t) {
    for (int i = 0; i != nclosed_*2; ++i) {
      const int nt = t + nclosed_*2; 
      (*mitti)(i, t) += (*cfock)(nt, nt) + (*afock)(nt, nt) - (*cfock)(i, i) - (*afock)(i, i)
                        + (*rdm1)(t,t)*(*cfock)(i, i) - (*cfockd)(t, t) - (*qxr->get_conjg())(nt ,t)
                        + (*miitt2rdm)(i, t) + (*miitt1rdm)(i, t);
    }                   
  } 

  if (verbose) {
    mitti->get_submatrix(morbital, norbital, 1, 1)->print("G(1,1)_{ti,ti}");
    cout << ">>>>>>>>>>>> debug : G(1,1)(ti,ti) >>>>>>>>>>>>" << endl << endl;
  }

  return mitti;
}


shared_ptr<ZMatrix> ZCASSCF::___debug___closed_active_diagonal_hessian_kramers(shared_ptr<const ZMatrix> coeffi, shared_ptr<const ZMatrix> coefft, const bool verbose) const {
  /* returns Mat(i,t) = G^{(1,1)}_{ti,ti} = [ (i ki|kt u) - (iu|kt i) ] D(t u) - [ (i ki|u t) - (u ki|i t) ] D(u kt) // CHECK SIGNS
                                          + [ (i ki|v u) - (i u|v ki) ] G(vu,t kt)
                                          + [ (kt ki|i t) - (i ki|kt t) ]
  for the time being, we implement it in the worst possible way... to be updated to make it efficient. */
  assert(coefft->mdim() == nact_*2);
  if (verbose)
    cout << ">>>>>>>>>>>> debug : G(1,1)(kt ki,t i) >>>>>>>>>>>>" << endl;

  shared_ptr<ZMatrix> kmitti1rdm  = ___debug___closed_active_diagonal_1rdm_contraction_exchange(coeffi, coefft); // (i t|u ki) * D(u kt)
  shared_ptr<ZMatrix> kmiitt1rdm  = ___debug___diagonal_1rdm_contraction_coulomb(coeffi, coefft, true); // (i ki| kt u) * D(tu) ; apears to be 0
  shared_ptr<ZMatrix> kmitti1rdmb = ___debug___diagonal_1rdm_contraction_exchange(coeffi, coefft, true); // (i u|kt ki) * D(t u) ; appears to be 0
  if(verbose) {
    kmitti1rdm->get_submatrix(morbital, norbital, 1, 1)->print("(u ki|i t) * D(u kt)");
    kmitti1rdmb->get_submatrix(morbital, norbital, 1, 1)->print("(i t|u ki) * D(u kt)");
    kmiitt1rdm->get_submatrix(morbital, norbital, 1, 1)->print("(i ki|kt u) * D(tu)");
  }

  shared_ptr<ZMatrix> kmitti = ___debug___diagonal_integrals_exchange_kramers(coeffi, coefft, false, true); // (i t|kt ki)
  shared_ptr<ZMatrix> kmiitt = ___debug___diagonal_integrals_coulomb_kramers(coeffi, coefft);  // (i ki|kt t) ; appears to be 0
  if (verbose) {
   kmitti->get_submatrix(morbital, norbital, 1, 1)->print("(i t|kt ki)");
   kmiitt->get_submatrix(morbital, norbital, 1, 1)->print("(i ki|kt t)");
  }
  *kmitti -= *kmiitt;
  *kmitti += (*kmiitt1rdm - *kmitti1rdmb - *kmitti1rdm); // TODO : still short one coulomb term, but it is probably 0 DOUBLE CHECK!

  shared_ptr<ZMatrix> krdm2coulomb = ___debug___diagonal_integrals_coulomb_active_kramers(coeffi, coefft, true); // appears 0
  shared_ptr<ZMatrix> krdm2exch    = ___debug___diagonal_2rdm_contraction_exchange(coeffi, true);
  if (verbose) {
    krdm2exch->get_submatrix(morbital, norbital, 1, 1)->print("(i u|v ki)*G(v u,t kt)");
    krdm2coulomb->get_submatrix(morbital, norbital, 1, 1)->print("(i ki|v u)*G(v u,t kt)");
  }
  *kmitti += (*krdm2coulomb - *krdm2exch);

  if (verbose) {
    kmitti->get_submatrix(morbital, norbital, 1, 1)->print("G(1,1)_{kt ki,t i}");
    cout << "<<<<<<<<<<<< debug : G(1,1)(kt ki,t i) <<<<<<<<<<<<" << endl << endl;
  }

  return kmitti;
}


shared_ptr<ZMatrix> ZCASSCF::___debug___closed_active_offdiagonal_hessian_kramers(shared_ptr<const ZMatrix> coeffi, shared_ptr<const ZMatrix> coefft, const bool verbose) const {
  /* returns Mat(i,t) = G^{(1,2)}_{ti,ti} = [ (t ki|v k) - (t i|v ki) ] D(v kt) + [ (kt i|v ki) - (kt ki|v i) ] D(v t)
                                          -   (i ki|v u)  G(t kt,v u)   + [ (t i|kt ki) - (t ki|kt i) ]
  for the time being, we implement it in the worst possible way... to be updated to make it efficient. */
  assert(coefft->mdim() == nact_*2);

  if (verbose)
    cout << "<<<<<<<<<<<< debug : G(1,2)(kt i,kt i) <<<<<<<<<<<<" << endl;

  shared_ptr<const ZMatrix> rdm1 = transform_rdm1();
  if (verbose)
    rdm1->print("1rdm");
  shared_ptr<ZMatrix> kmitit     = ___debug___diagonal_integrals_exchange_kramers(coeffi, coefft, true); 
  shared_ptr<ZMatrix> offd1rdmx  = ___debug___closed_active_offdiagonal_1rdm_exchange(coeffi, coefft);
  shared_ptr<ZMatrix> kmitit2rdm = ___debug___closed_active_offdiagonal_2rdm_exchange(coeffi, coefft);
  if (verbose) {
    kmitit->get_submatrix(morbital, norbital, 1, 1)->print("(t i|kt ki) - (t ki|kt i)");
    kmitit2rdm->get_submatrix(morbital, norbital, 1, 1)->print("(v i| u ki) * G(u t, v kt)");
    offd1rdmx->get_submatrix(morbital, norbital, 1, 1)->print("(t ki|v i) D(v kt) + ...");
  }
  *kmitit += (*offd1rdmx - *kmitit2rdm);

  if (verbose) {
    kmitit->get_submatrix(morbital, norbital, 1, 1)->print("G(1,2)_{kt i,kt i}");
    cout << ">>>>>>>>>>>> debug : G(1,2)(kt i,kt i) >>>>>>>>>>>>" << endl << endl;
  }

  return kmitit;
}
