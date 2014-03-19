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
  // todo can remove when have task
  deninact_ = cas->ao_rdm1(fci_->rdm1(0), true);
  denall_ = cas->ao_rdm1(fci_->rdm1(0));

  // for denominator
  eig_ = assemble_eig(ref_->coeff(), ref_->hcore());

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
#if 0 // TODO check
  shared_ptr<const Coeff> coeff = smith->coeff();
#else
  shared_ptr<const Coeff> coeff = ref_->coeff();
#endif

  // compute Yrs
  compute_y(d1, correction, d2, cider, coeff);
  yrs_->print("Yrs in IN GRAD mo basis", 20);

  // solve CPCASSCF
  auto g0 = make_shared<Matrix>(*yrs_);
  auto g1 = make_shared<Dvec>(cider,1);
  auto grad = make_shared<PairFile<Matrix, Dvec>>(g0, g1);

  shared_ptr<DFHalfDist> half   = ref_->geom()->df()->compute_half_transform(coeff->slice(0, ref_->nocc()));
  shared_ptr<DFHalfDist> halfj  = half->apply_J();
  shared_ptr<DFHalfDist> halfjj = halfj->apply_J();


#if 0
  assert(check_blocks(g0));
  auto cp = make_shared<CPCASSCF>(grad, fci_->civectors(), half, halfjj, ref_, fci_);
  shared_ptr<const Matrix> zmat, xmat;
  shared_ptr<const Dvec> zvec;
  tie(zmat, zvec, xmat) = cp->solve();

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

shared_ptr<Matrix> CASPT2Grad::assemble_eig(shared_ptr<const Coeff> coeff, shared_ptr<const Matrix> hcore) {
  const int target = 0;
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();
  const int nvirt = ref_->nvirt();

  // related to denominators
  const int nmobasis = coeff->mdim();
  assert(nmobasis == nclosed+nact+nvirt);
  auto eig = make_shared<Matrix>(nmobasis, nmobasis);
  {
    // as in Theor Chem Acc (1997) 97:88-95
    vector<double> occup_ = fci_->rdm1(target)->diag();

    auto finact = make_shared<Matrix>(*coeff % *(fci_->jop()->core_fock()) * *coeff);

    auto denact = make_shared<Matrix>(*denall_-*deninact_);
    auto fact_ao = make_shared<Fock<1>>(geom_, hcore, denact, ref_->schwarz());
    auto f = make_shared<Matrix>(*finact+ *coeff%(*fact_ao-*hcore)**coeff);


#if 1 // TODO determine proper rdm
    auto fact = make_shared<Qvec>(nmobasis, nact, coeff, nclosed, fci_, fci_->rdm2(target));
#else
    auto fact = make_shared<Qvec>(nmobasis, nact, coeff, nclosed, fci_, fci_->rdm2_av());
#endif
    for (int i = 0; i != nact; ++i)
      daxpy_(nmobasis, occup_[i], finact->element_ptr(0,nclosed+i), 1, fact->data()+i*nmobasis, 1);

    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nvirt; ++j)
        eig->element(j+nocc,i+nclosed) = eig->element(i+nclosed,j+nocc) = -fact->element(i,i) + occup_[i]*f->element(j+nocc, j+nocc);

    for (int i = 0; i != nclosed; ++i)
      for (int j = 0; j != nvirt; ++j)
         eig->element(j+nocc,i) = eig->element(i,j+nocc) = 2.0*f->element(j+nocc, j+nocc) - 2.0*f->element(i, i);

    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nclosed; ++j) {
         eig->element(j,i+nclosed) = eig->element(i+nclosed,j)
#if 1 // original
                                   = (f->element(nclosed+i,nclosed+i)*2.0-fact->element(i+nclosed,i)) - f->element(j, j)*(2.0 - occup_[i]);
#else // this converges cpcasscf A11
                                   = (f->element(nclosed+i,nclosed+i)*1.0-fact->element(i+nclosed,i)) - f->element(j, j)*(2.0 - occup_[i]);
#endif
      }

  }

  return eig;
}


void CASPT2Grad::compute_y(shared_ptr<const Matrix> dm1, double correction, shared_ptr<const Matrix> dm2, shared_ptr<const Civec> cider, shared_ptr<const Coeff> coeff_) {
  const int target = 0;
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
  // add correction to active space
  for (int i = nclosed; i != nclosed+nact; ++i) dmr->element(i, i) -=  correction * 2.0;

  // TODO they are redundant, though...
  shared_ptr<DFHalfDist> half   = ref_->geom()->df()->compute_half_transform(coeff_->slice(0, nocc));
  shared_ptr<DFHalfDist> halfj  = half->apply_J();
  shared_ptr<DFHalfDist> halfjj = halfj->apply_J();
  shared_ptr<const DFFullDist> full = halfj->compute_second_transform(coeff_);

  // Y_rs = 2[Y1 + Y2 + Y3(ri) + Y4 + Y5(ri)]
  shared_ptr<Matrix> out = make_shared<Matrix>(nmobasis, nmobasis);
  auto y    = make_shared<Matrix>(nmobasis, nmobasis);
  auto y1   = make_shared<Matrix>(nmobasis, nmobasis);
  auto y2   = make_shared<Matrix>(nmobasis, nmobasis);
  auto y3ri = make_shared<Matrix>(nmobasis, nocc);
  auto y4   = make_shared<Matrix>(nmobasis, nmobasis);
  auto y5ri = make_shared<Matrix>(nmobasis, nocc);

  {
    // 2 Y1 = h(d0 + d2) * 2
    // one-electron contributions
    auto hmo = make_shared<const Matrix>(*coeff_ % *ref_->hcore() * *coeff_);
    auto dtot = make_shared<Matrix>(*dmr);
    for (int i = 0; i != nclosed; ++i) dtot->element(i, i) += 2.0;
    dtot->add_block(1.0, nclosed, nclosed, nact, nact, ref_->rdm1(target)->data());
    *y1 = *hmo * *dtot * 2.0;
  }

  {
    // Y2 = Y2_rs = Y2_ri + Y2_ra, so making both at once
    // 2 Y2_rs = 2[[(rt|kl) -1/2(rk|tl)]d^(0)_{kl}]*dm1_ts
    shared_ptr<const Matrix> dkl = ref_->rdm1_mat(target);
    auto ocoeff = coeff_->slice(0, nocc);
    auto dklao = make_shared<const Matrix>(*ocoeff * *dkl ^ *ocoeff);
    // J_rt = (rt|delta) (gamma|kl)*d_kl
    auto jrt = make_shared<const Matrix>(*coeff_ % *(ref_->geom())->df()->compute_Jop(dklao) * *coeff_);
    // 2 Y2_rs += 2*J_rt*d2_ts
    *y2 = *jrt * *dmr * 2.0;
    //  -1 K_tr  = -(rk|lt)d_kl
    // make (G'|kt)
    auto halfj_kt = halfj->copy();
    halfj_kt->rotate_occ(dkl);
    shared_ptr<const Matrix> mat = halfj->form_2index(halfj_kt, -1.0);
    auto krt = make_shared<const Matrix>(*coeff_ % *mat * *coeff_);
    // Y2_rs += Krt*d2_ts
    *y2 += *krt * *dmr;
  }

  {
    // 2 Y3 = 2 Y3_ri = 2[sum_st [(rj|st) -1/2(rs|tj)]d^(2)_{st}]*dm0_ji, where st can be cc,xx,aa
    // convert d2 to ao basis
    auto dmrao = make_shared<const Matrix>(*coeff_ * *dmr ^ *coeff_);
    dmrao->print("printing dmr2ao caspt2", 20);
    //  2 Jrj = 2 (rj|st)d2_st
    auto jrj = make_shared<const Matrix>(*coeff_ % *(ref_->geom())->df()->compute_Jop(dmrao) * *ocmat);
    *y3ri = *jrj * *(ref_->rdm1_mat(0)) * 2.0;
    // -1 K_jr  = -(rs|jt)d2_st
    auto kjr = make_shared<const Matrix>(*(halfjj->compute_Kop_1occ(dmrao, -1.0)) * *coeff_);
    *y3ri += *kjr % *(ref_->rdm1_mat(0));
  }

  {
    // make 2 Y4 = 2 y4_rs = 2[y4ri + y4ra] = 2 K^{kl}_{rt} D^{lk}_{ts} = 2 (kr|lt) D_(lt,ks) = 2 (kr|lt)[D0_(lt,ks) + D1_(lt,ks)]
    auto Dtot = make_shared<Matrix>(nocc*nall, nocc*nall);
    auto D1 = make_shared<Matrix>(nocc*nall, nocc*nall);
    // extend 2rdm dimensions from D0(lj,ki) -> D0(lt,ks)
    shared_ptr<const Matrix> D0 = ref_->rdm2_mat_2gen(0);
    // resizing dm2_(le,kf) to dm2_(ls,kt) and saving as D1_(lt,ks)
    for (int s = 0; s != nall; ++s) // extend
      for (int k = 0; k != nocc; ++k)
        for (int t = 0; t != nall; ++t) // extend
          for (int l = 0; l != nocc; ++l) {
            if (t >= nclosed && s >= nclosed) {
              D1->element(l+nocc*s, k+nocc*t) = dm2->element(l+nocc*(t-nclosed), k+nocc*(s-nclosed));
            }
          }

    *Dtot+=*D0;
    *Dtot+=*D1;
    // form (G'|ks)
    // todo check transpose
    shared_ptr<const DFFullDist> fullks = full->apply_2rdm(Dtot->data());
    y4 = full->form_2index(fullks, 2.0);
  }

  {
    // 2 Y5 = 2 Y5_ri = 2 Ybar (Gyorffy)  = 2 (rs|tj) D^ij_st = 2 (rs|jt) [D0(is,jt) + D1(is,jt)]
    // D0 needs to be extended from D0(il,jk) to D0(is,jt)
    shared_ptr<const Matrix> D0st = ref_->rdm2_mat_2gen(0);

    // need to make Dtot_st
    auto Dtot_st = make_shared<Matrix>(nocc*nall, nocc*nall);
    auto D1st = make_shared<Matrix>(nocc*nall, nocc*nall);
    // resize dm2_(le,kf) to lt,ks then save as ls,kt
    for (int i = 0; i != nall; ++i) // extending
      for (int k = 0; k != nocc; ++k)
        for (int s = 0; s != nall; ++s) // extending
          for (int l = 0; l != nocc; ++l) {
            if (s >= nclosed &&  i >= nclosed) {
              D1st->element(l+nocc*i, k+nocc*s) = dm2->element(l+nocc*(s-nclosed), k+nocc*(i-nclosed));
            }
          }

    *Dtot_st += *D0st;
    *Dtot_st += *D1st;

    // todo check transpose
    shared_ptr<const DFFullDist> fullis = full->apply_2rdm(Dtot_st->data());
    shared_ptr<const DFHalfDist> dfback = fullis->back_transform(coeff_)->apply_J();
    auto y5ri_ao = ref_->geom()->df()->form_2index(dfback, 2.0);
    *y5ri = *coeff_ % *y5ri_ao;
  }

  // make Yrs in mo basis
  *y += *y1;
  *y += *y2;
  y->add_block(1.0, 0, 0, nmobasis, nocc, y3ri);
  *y += *y4;
  y->add_block(1.0, 0, 0, nmobasis, nocc, y5ri);

  yrs_ = y;
}


