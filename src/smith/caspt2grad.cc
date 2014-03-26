//
// BAGEL - Parallel electron correlation program.
// Filename: caspt2grad.cc
// Copyright (C) 2013 Toru Shiozaki
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


#include <src/smith/caspt2grad.h>
#include <src/casscf/superci.h>
#include <src/casscf/qvec.h>
#include <src/smith/smith.h>
#include <src/grad/gradeval_base.h>
#include <src/math/algo.h>


using namespace std;
using namespace bagel;

CASPT2Grad::CASPT2Grad(shared_ptr<const PTree> inp, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
  : Method(inp, geom, ref) {

  // compute CASSCF first
  auto cas = make_shared<SuperCI>(inp, geom, ref);
  cas->compute();

  cout << endl << "  === DF-CASPT2Grad calculation ===" << endl << endl;
  if (geom->df() == nullptr) throw logic_error("CASPT2Grad is only implemented with DF");

  // update reference
  ref_ = cas->conv_to_ref();
  fci_ = cas->fci();
}


void CASPT2Grad::compute() {
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();

  // state-averaged density matrices
  shared_ptr<const RDM<1>> rdm1_av = fci_->rdm1_av();
  shared_ptr<const RDM<2>> rdm2_av = fci_->rdm2_av();

  shared_ptr<const Matrix> d1, d2;
  shared_ptr<const Civec> cider;
  {
    // construct SMITH here
    shared_ptr<const PTree> smithinput = idata_->get_child("smith");
    auto smith = make_shared<Smith>(smithinput, ref_->geom(), ref_);
    smith->compute();

    // use coefficients from smith (closed and virtual parts have been rotated in smith to make them canonical).
    coeff_ = smith->coeff();

    // save correlated density matrices d(1), d(2), and ci derivatives
    shared_ptr<Matrix> d1tmp = make_shared<Matrix>(*smith->dm1());
    const double correction = smith->correction();
    // add correction to active part of the correlated one-body density
    for (int i = nclosed; i != nclosed+nact; ++i)
      d1tmp->element(i, i) -=  correction * 2.0;
    d1 = d1tmp;
    d2 = smith->dm2();
    cider = smith->cider();
    target_ = smith->algo()->ref()->target();
  }

  const int nmobasis = coeff_->mdim();
  shared_ptr<const Matrix> d0 = ref_->rdm1_mat(target_)->resize(nmobasis,nmobasis);
  shared_ptr<const Matrix> ocoeff = coeff_->slice(0, nocc);

  {
    auto dtotao = make_shared<Matrix>(*coeff_ * (*d0 + *d1) ^ *coeff_);
    Dipole dipole(geom_, dtotao, "CASPT2 unrelaxed");
  }

  // compute Yrs
  shared_ptr<const DFHalfDist> half   = ref_->geom()->df()->compute_half_transform(ocoeff);
  shared_ptr<const DFHalfDist> halfj  = half->apply_J();
  shared_ptr<const DFHalfDist> halfjj = halfj->apply_J();
  shared_ptr<Matrix> yrs;
  shared_ptr<const DFFullDist> fulld1; // (gamma| ir) D(ir,js)
  tie(yrs, fulld1) = compute_y(d1, d2, cider, half, halfj, halfjj);

  // solve CPCASSCF
  auto g0 = yrs;
  auto g1 = make_shared<Dvec>(cider, ref_->nstate()); // FIXME this is wrong for nstate > 1 in CASSCF
  auto grad = make_shared<PairFile<Matrix, Dvec>>(g0, g1);

  auto cp = make_shared<CPCASSCF>(grad, fci_->civectors(), half, halfjj, ref_, fci_, coeff_);
  shared_ptr<const Matrix> zmat, xmat;
  shared_ptr<const Dvec> zvec;
  tie(zmat, zvec, xmat) = cp->solve(1.0e-10);

  // form relaxed 1RDM
  // form Zd + dZ^+
  shared_ptr<Matrix> dsa = rdm1_av->rdm1_mat(nclosed)->resize(nmobasis, nmobasis);
  auto dm = make_shared<Matrix>(*zmat * *dsa + (*dsa ^ *zmat));

  shared_ptr<Matrix> dtot = d0->copy();
  dtot->ax_plus_y(1.0, dm);
  dtot->ax_plus_y(1.0, d1);

  // form zdensity
  auto detex = make_shared<Determinants>(fci_->norb(), fci_->nelea(), fci_->neleb(), false, /*mute=*/true);
  shared_ptr<const RDM<1>> zrdm1;
  shared_ptr<const RDM<2>> zrdm2;
  tie(zrdm1, zrdm2) = fci_->compute_rdm12_av_from_dvec(fci_->civectors(), zvec, detex);

  shared_ptr<Matrix> zrdm1_mat = zrdm1->rdm1_mat(nclosed, false)->resize(nmobasis, nmobasis);
  zrdm1_mat->symmetrize();
  dtot->ax_plus_y(1.0, zrdm1_mat);

  // compute relaxed dipole to check
  auto dtotao = make_shared<Matrix>(*coeff_ * *dtot ^ *coeff_);
  {
    Dipole dipole(geom_, dtotao, "CASPT2 relaxed");
    dipole.compute();
  }

  // xmat in the AO basis
  auto xmatao = make_shared<Matrix>(*coeff_ * *xmat ^ *coeff_);

  // two-body part
  // first make occ-occ part (copy-and-paste from src/casscf/supercigrad.cc)
  shared_ptr<const DFFullDist> qij  = halfjj->compute_second_transform(ocoeff);
  shared_ptr<DFHalfDist> qri;
  {
    shared_ptr<const Matrix> ztrans = make_shared<Matrix>(*coeff_ * *zmat->slice(0,nocc));
    {
      const RDM<2> D(*ref_->rdm2(target_)+*zrdm2);
      const RDM<1> dd(*ref_->rdm1(target_)+*zrdm1);

      shared_ptr<DFFullDist> qijd = qij->apply_2rdm(D.data(), dd.data(), nclosed, nact);
      qijd->ax_plus_y(2.0, halfjj->compute_second_transform(ztrans)->apply_2rdm(rdm2_av->data(), rdm1_av->data(), nclosed, nact));
      qri = qijd->back_transform(ocoeff);
    }
    {
      shared_ptr<const DFFullDist> qijd2 = qij->apply_2rdm(rdm2_av->data(), rdm1_av->data(), nclosed, nact);
      qri->ax_plus_y(2.0, qijd2->back_transform(ztrans));
    }
  }

  // D1 part. 2.0 seems to come from the difference between smith and bagel (?)
  qri->ax_plus_y(2.0, fulld1->apply_J()->back_transform(coeff_));

  // contributions from non-separable part
  shared_ptr<Matrix> qq  = qri->form_aux_2index(halfjj, 1.0);
  shared_ptr<DFDist> qrs = qri->back_transform(ocoeff);

  // separable part
  // size of naux
  {
    shared_ptr<const Matrix> d0ao = make_shared<Matrix>(*coeff_ * *d0 ^ *coeff_);
    shared_ptr<const Matrix> d1ao = make_shared<Matrix>(*coeff_ * *d1 ^ *coeff_);
    shared_ptr<const Matrix> cd0 = geom_->df()->compute_cd(d0ao);
    shared_ptr<const Matrix> cd1 = geom_->df()->compute_cd(d1ao);

    // three-index derivatives (seperable part)...
    vector<shared_ptr<const Matrix>> cd {cd0, cd1};
    vector<shared_ptr<const Matrix>> dd {d1ao, d0ao};

    shared_ptr<DFHalfDist> sepd = halfjj->apply_density(d1ao);
    sepd->rotate_occ(ref_->rdm1_mat(target_)); // d0 in occ-occ size

    qrs->ax_plus_y(-1.0, sepd->back_transform(ocoeff)); // TODO this back transformation can be done together
    qrs->add_direct_product(cd, dd, 1.0);

    *qq += (*cd0 ^ *cd1) * 2.0;
    *qq += *halfjj->form_aux_2index(sepd, -1.0);
  }

  // compute gradients
  GradEval_base g(geom_);
  shared_ptr<GradFile> gradient = g.contract_gradient(dtotao, xmatao, qrs, qq);
  gradient->print();
}


tuple<shared_ptr<Matrix>, shared_ptr<const DFFullDist>>
  CASPT2Grad::compute_y(shared_ptr<const Matrix> dm1, shared_ptr<const Matrix> dm2, shared_ptr<const Civec> cider,
                        shared_ptr<const DFHalfDist> half, shared_ptr<const DFHalfDist> halfj, shared_ptr<const DFHalfDist> halfjj) {
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();
  const int nvirt = ref_->nvirt();
  const int nall = nocc + nvirt;
  const int nmobasis = coeff_->mdim();
  assert(nall == nmobasis);

  shared_ptr<const Matrix> ocmat = coeff_->slice(0, nocc);
  shared_ptr<const Matrix> vcmat = coeff_->slice(nocc, nmobasis);

  auto dmr = make_shared<Matrix>(*dm1);

  shared_ptr<const DFFullDist> full = halfj->compute_second_transform(coeff_);
  shared_ptr<const DFFullDist> fullo = halfj->compute_second_transform(ocmat);

  // Y_rs = 2[Y1 + Y2 + Y3(ri) + Y4 + Y5(ri)]
  shared_ptr<Matrix> out = make_shared<Matrix>(nmobasis, nmobasis);
  auto ocoeff = coeff_->slice(0, nocc);

  {
    // 2 Y1 = h(d0 + d1 + d2) * 2
    // one-electron contributions
    auto hmo = make_shared<const Matrix>(*coeff_ % *ref_->hcore() * *coeff_);
    auto d0 = make_shared<Matrix>(*ref_->rdm1_mat(target_)->resize(nmobasis,nmobasis));
    *out += *hmo * (*dmr + *d0) * 2.0;
  }

  {
    // Y2 = Y2_rs = Y2_ri + Y2_ra, so making both at once
    shared_ptr<Matrix> dkl = ref_->rdm1_mat(target_);
    dkl->sqrt();
    dkl->scale(1.0/sqrt(2.0));
    Fock<1> fock(geom_, ref_->hcore()->clone(), nullptr, make_shared<Matrix>(*ocoeff * *dkl), /*grad*/false, /*rhf*/true);
    *out += *coeff_ % fock * *coeff_ * *dmr * 2.0;
  }

  {
    // 2 Y3 = 2 Y3_ri*dm0_ji
    // coulomb
    auto dmrao = make_shared<Matrix>(*coeff_ * *dmr ^ *coeff_);
    auto jop = geom_->df()->compute_Jop(dmrao);
    // exchange
    auto kopi = halfjj->compute_Kop_1occ(dmrao, -0.5)->transpose();

    out->add_block(2.0, 0, 0, nmobasis, nocc, *coeff_ % (*jop * *ocoeff + *kopi) * *ref_->rdm1_mat(target_));
  }

  // TODO D1 must be parallelised as it is very big.
  // construct D1 to be used in Y4 and Y5
  auto D1 = make_shared<Matrix>(nocc*nall, nocc*nall);
  {
    // resizing dm2_(le,kf) to dm2_(lt,ks). no resort necessary.
    for (int s = 0; s != nall; ++s) // extend
      for (int k = 0; k != nocc; ++k)
        for (int t = 0; t != nall; ++t) // extend
          for (int l = 0; l != nocc; ++l) {
            if (t >= nclosed && s >= nclosed) {
              D1->element(l+nocc*t, k+nocc*s) = dm2->element(l+nocc*(t-nclosed), k+nocc*(s-nclosed));
            }
          }
  }

  // fullks will be reused for gradient contraction
  shared_ptr<const DFFullDist> fullks;
  {
    // 2 Y4 =  2 K^{kl}_{rt} D^{lk}_{ts} = 2 (kr|lj) D0_(lj,ki) +  2 (kr|lt) D1_(lt,ks)
    // construct stepwise, D1 part
    fullks = full->apply_2rdm(D1->data());
    *out += *full->form_2index(fullks, 2.0);
    // D0 part
    shared_ptr<const DFFullDist> fulld = fullo->apply_2rdm(ref_->rdm2(target_)->data(), ref_->rdm1(target_)->data(), nclosed, nact);
    out->add_block(2.0, 0, 0, nmobasis, nocc, full->form_2index(fulld, 1.0));
  }

  {
    // 2 Y5 = 2 Y5_ri = 2 Ybar (Gyorffy)  = 2 (rs|tj) D^ij_st = 2 (rl|jk) D0_(il,jk) + 2 (rs|tj) D1_(is,jt)]
    // construct stepwise, D1 part
    shared_ptr<const DFHalfDist> dfback = fullks->apply_J()->back_transform(coeff_);
    auto y5ri_ao = ref_->geom()->df()->form_2index(dfback, 1.0);
    out->add_block(2.0, 0, 0, nmobasis, nocc, *coeff_ % *y5ri_ao);
  }

  return make_tuple(out, fullks);
}

