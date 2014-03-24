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


using namespace std;
using namespace bagel;

CASPT2Grad::CASPT2Grad(shared_ptr<const PTree> inp, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref) : Method(inp, geom, ref),  target_state_(inp->get<int>("target", 0)) {

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
  // construct SMITH here
  shared_ptr<const PTree> smithinput = idata_->get_child("smith");
  auto smith = make_shared<Smith>(smithinput, ref_->geom(), ref_);

  smith->compute();

  // save correlated density matrices d(1), d(2), and ci derivatives
  shared_ptr<const Matrix> d1 = smith->dm1();
  shared_ptr<const Matrix> d2 = smith->dm2();
  double correction = smith->correction();
  shared_ptr<const Civec> cider = smith->cider();
  shared_ptr<const Coeff> coeff = smith->coeff();

  {
    const int nmobasis = coeff->mdim();
    auto dtotao = make_shared<Matrix>(*coeff * *ref_->rdm1_mat(target_state_)->resize(nmobasis,nmobasis) ^ *coeff);
    Dipole dipole(geom_, dtotao, "CASPT2 unrelaxed");
  }

  // compute Yrs
  auto yrs = compute_y(d1, correction, d2, cider, coeff);
  yrs->antisymmetrize();
  yrs->print("Yrs in IN GRAD mo basis");

  // solve CPCASSCF
  auto g0 = yrs;
auto tmp = yrs->copy();
tmp->antisymmetrize();
tmp->print();
  auto g1 = make_shared<Dvec>(cider, ref_->nstate()); // FIXME this is wrong for nstate > 1 in CASSCF
  auto grad = make_shared<PairFile<Matrix, Dvec>>(g0, g1);

  shared_ptr<DFHalfDist> half   = ref_->geom()->df()->compute_half_transform(coeff->slice(0, ref_->nocc()));
  shared_ptr<DFHalfDist> halfjj = half->apply_JJ();

#if 1
  //assert(check_blocks(g0));
  auto cp = make_shared<CPCASSCF>(grad, fci_->civectors(), half, halfjj, ref_, fci_, coeff);
  shared_ptr<const Matrix> zmat, xmat;
  shared_ptr<const Dvec> zvec;
  tie(zmat, zvec, xmat) = cp->solve(1.0e-10);

  // form relaxed 1RDM
  // form Zd + dZ^+
  const int target = target_state_;
  const int nclosed = ref_->nclosed();
  const int nmobasis = coeff->mdim();
  shared_ptr<Matrix> dsa = fci_->rdm1_av()->rdm1_mat(nclosed)->resize(nmobasis, nmobasis);
  auto dm = make_shared<Matrix>(*zmat * *dsa + (*dsa ^ *zmat));

  shared_ptr<Matrix> dtot = ref_->rdm1_mat(target)->resize(nmobasis, nmobasis);
  dtot->ax_plus_y(1.0, dm);

  // form zdensity
  auto detex = make_shared<Determinants>(fci_->norb(), fci_->nelea(), fci_->neleb(), false, /*mute=*/true);
  shared_ptr<const RDM<1>> zrdm1;
  shared_ptr<const RDM<2>> zrdm2;
  tie(zrdm1, zrdm2) = fci_->compute_rdm12_av_from_dvec(fci_->civectors(), zvec, detex);

  shared_ptr<Matrix> zrdm1_mat = zrdm1->rdm1_mat(nclosed, false)->resize(nmobasis, nmobasis);
  zrdm1_mat->symmetrize();
  dtot->ax_plus_y(1.0, zrdm1_mat);

  // compute relaxed dipole to check
  auto dtotao = make_shared<Matrix>(*coeff * *dtot ^ *coeff);
  Dipole dipole(geom_, dtotao, "CASPT2 relaxed");
  dipole.compute();
#endif

  // form relaxed 2RDM

  // compute gradients

}

bool CASPT2Grad::check_blocks(shared_ptr<Matrix> grad) {
  // form Y =  Yrs - Ysr
  grad->antisymmetrize();
  assert(grad->ndim() == grad->mdim());
  bool out = true;
  // closed-closed block should be zero
  for (int j = 0;  j != ref_->nclosed(); ++j) {
    for (int i = 0;  i != ref_->nclosed(); ++i) {
      if (fabs(grad->element(i,j)) > numerical_zero__) {
        out = false;
        break;
      }
    }
    if (out == false) break;
  }

  // virtual-virtual block should be zero
  if (out != false)
  for (int j = 0;  j != ref_->nvirt(); ++j) {
    for (int i = 0;  i != ref_->nvirt(); ++i) {
      if (fabs(grad->element(i+ref_->nocc(),j+ref_->nocc())) > numerical_zero__) {
        out = false;
        break;
      }
    }
    if (out == false) break;
  }

  return out;
}


shared_ptr<Matrix> CASPT2Grad::compute_y(shared_ptr<const Matrix> dm1, double correction, shared_ptr<const Matrix> dm2, shared_ptr<const Civec> cider, shared_ptr<const Coeff> coeff_) {
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
  // add correction to active part of the correlated one-body density
  for (int i = nclosed; i != nclosed+nact; ++i) dmr->element(i, i) -=  correction * 2.0;

  // TODO they are redundant, though...
  shared_ptr<DFHalfDist> half   = ref_->geom()->df()->compute_half_transform(ocmat);
  shared_ptr<DFHalfDist> halfj  = half->apply_J();
  shared_ptr<DFHalfDist> halfjj = halfj->apply_J();
  shared_ptr<const DFFullDist> full = halfj->compute_second_transform(coeff_);
  shared_ptr<const DFFullDist> fullo = halfj->compute_second_transform(ocmat);

  // Y_rs = 2[Y1 + Y2 + Y3(ri) + Y4 + Y5(ri)]
  shared_ptr<Matrix> out = make_shared<Matrix>(nmobasis, nmobasis);
  auto ocoeff = coeff_->slice(0, nocc);

  {
    // 2 Y1 = h(d0 + d1 + d2) * 2
    // one-electron contributions
    auto hmo = make_shared<const Matrix>(*coeff_ % *ref_->hcore() * *coeff_);
    auto d0 = make_shared<Matrix>(*ref_->rdm1_mat(target_state_)->resize(nmobasis,nmobasis));
    *out += *hmo * (*dmr + *d0) * 2.0;
  }

  {
    // Y2 = Y2_rs = Y2_ri + Y2_ra, so making both at once
    shared_ptr<Matrix> dkl = ref_->rdm1_mat(target_state_);
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

    out->add_block(2.0, 0, 0, nmobasis, nocc, *coeff_ % (*jop * *ocoeff + *kopi) * *ref_->rdm1_mat(target_state_));
  }

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

  {
    // 2 Y4 =  2 K^{kl}_{rt} D^{lk}_{ts} = 2 (kr|lj) D0_(lj,ki) +  2 (kr|lt) D1_(lt,ks)
    // construct stepwise, D1 part
    shared_ptr<const DFFullDist> fullks = full->apply_2rdm(D1->data());
    *out += *full->form_2index(fullks, 2.0);
    // D0 part
    shared_ptr<const DFFullDist> fulld = fullo->apply_2rdm(ref_->rdm2(target_state_)->data(), ref_->rdm1(target_state_)->data(), nclosed, nact);
    out->add_block(2.0, 0, 0, nmobasis, nocc, full->form_2index(fulld, 1.0));
  }

  {
    // 2 Y5 = 2 Y5_ri = 2 Ybar (Gyorffy)  = 2 (rs|tj) D^ij_st = 2 (rl|jk) D0_(il,jk) + 2 (rs|tj) D1_(is,jt)]
    // construct stepwise, D1 part
    shared_ptr<const DFFullDist> fullis = full->apply_2rdm(D1->data());
    shared_ptr<const DFHalfDist> dfback = fullis->back_transform(coeff_)->apply_J();
    auto y5ri_ao = ref_->geom()->df()->form_2index(dfback, 1.0);
    out->add_block(2.0, 0, 0, nmobasis, nocc, *coeff_ % *y5ri_ao);
    // D0 part
    shared_ptr<const DFFullDist> fulljk = fullo->apply_2rdm(ref_->rdm2(target_state_)->data(), ref_->rdm1(target_state_)->data(), nclosed, nact);
    out->add_block(2.0, 0, 0, nmobasis, nocc, full->form_2index(fulljk, 1.0));
  }

  return out;
}


