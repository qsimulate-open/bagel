//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: caspt2grad.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bagel_config.h>

#include <src/scf/hf/fock.h>
#include <src/grad/cpcasscf.h>
#include <src/grad/gradeval.h>
#include <src/multi/casscf/cassecond.h>
#include <src/multi/casscf/casnoopt.h>
#include <src/multi/casscf/qvec.h>
#include <src/smith/smith.h>
#include <src/smith/caspt2grad.h>
#include <src/prop/multipole.h>
#include <src/prop/hyperfine.h>


using namespace std;
using namespace bagel;

CASPT2Grad::CASPT2Grad(shared_ptr<const PTree> inp, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
  : CASPT2Deriv(inp, geom, ref) {
#ifdef COMPILE_SMITH
  Timer timer;

  // compute CASSCF first
  if (inp->get<string>("algorithm", "") != "noopt") {
    auto cas = make_shared<CASSecond>(inp, geom, ref);
    cas->compute();
    ref_ = cas->conv_to_ref();
    fci_ = cas->fci();
    thresh_ = cas->thresh();
  } else {
    auto cas = make_shared<CASNoopt>(inp, geom, ref);
    cas->compute();
    ref_ = cas->conv_to_ref();
    fci_ = cas->fci();
    thresh_ = cas->thresh();
  }

  // gradient/property calculation
  target_ = inp->get<int>("_target");
  do_hyperfine_ = inp->get<bool>("hyperfine", false);

  maxziter_ = inp->get<int>("_maxziter");

  timer.tick_print("Reference calculation");

  cout << endl << "  === DF-CASPT2Grad calculation ===" << endl << endl;
#else
  throw logic_error("CASPT2 gradients require SMITH-generated code. Please compile BAGEL with --enable-smith");
#endif
}


// compute smith and set rdms and ci deriv to a member
void CASPT2Grad::compute() {
#ifdef COMPILE_SMITH
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();

  // construct SMITH here
  shared_ptr<PTree> smithinput = idata_->get_child("smith");
  smithinput->put("_grad", true);
  smithinput->put("_target", target_);
  smithinput->put("_hyperfine", do_hyperfine_);
  auto smith = make_shared<Smith>(smithinput, ref_->geom(), ref_);
  smith->compute();

  // use coefficients from smith (closed and virtual parts have been rotated in smith to make them canonical).
  coeff_ = smith->coeff();

  if (nact)
    cideriv_ = smith->cideriv()->copy();
  ncore_   = smith->algo()->info()->ncore();
  wf1norm_ = smith->wf1norm();
  msrot_   = smith->msrot();
  nstates_ = wf1norm_.size();
  assert(msrot_->ndim() == nstates_ && msrot_->mdim() == nstates_);
  assert(nstates_ == smith->algo()->info()->ciwfn()->nstates());

  Timer timer;

  // save correlated density matrices d(1), d(2), and ci derivatives
  auto d1tmp = make_shared<Matrix>(*smith->dm1());
  auto d11tmp = make_shared<Matrix>(*smith->dm11());
  d1tmp->symmetrize();
  d11tmp->symmetrize();
  // a factor of 2 from the Hylleraas functional (which is not included in the generated code)
  d11tmp->scale(2.0);

  auto d1set = [this](shared_ptr<const Matrix> d1t) {
    if (!ncore_) {
      return d1t->copy();
    } else {
      auto out = make_shared<Matrix>(coeff_->mdim(), coeff_->mdim());
      out->copy_block(ncore_, ncore_, coeff_->mdim()-ncore_, coeff_->mdim()-ncore_, d1t);
      return out;
    }
  };
  d1_ = d1set(d1tmp);
  d11_ = d1set(d11tmp);

  // XMS density matrix
  if (smith->dcheck()) {
    shared_ptr<const Matrix> dc = smith->dcheck();
    assert(dc->ndim() == nact && dc->mdim() == nact);
    auto tmp = make_shared<Matrix>(nocc, nocc);
    tmp->add_block(1.0, nclosed, nclosed, nact, nact, dc);
    dcheck_ = tmp;
  }

  // spin density matrices
  if (do_hyperfine_) {
    auto sd11tmp = make_shared<Matrix>(*smith->sdm11());
    sd11tmp->symmetrize();
    sd11tmp->scale(2.0);
    sd1_ = d1set(smith->sdm1());
    sd11_ = d1set(sd11tmp);
  }

  // zeroth order RDM
  if (nact) {
    d10ms_ = make_shared<RDM<1>>(nact);
    d20ms_ = make_shared<RDM<2>>(nact);
    for (int ist = 0; ist != nstates_; ++ist) {
      const double ims = msrot(ist, target());
      for (int jst = 0; jst != nstates_; ++jst) {
        const double jms = msrot(jst, target());
        shared_ptr<const RDM<1>> rdm1t;
        shared_ptr<const RDM<2>> rdm2t;
        tie(rdm1t, rdm2t) = ref_->rdm12(jst, ist, /*recompute*/true);
        d10ms_->ax_plus_y(ims*jms, *rdm1t);
        d20ms_->ax_plus_y(ims*jms, *rdm2t);
      }
    }
  }
  d2_ = smith->dm2();
  energy_ = smith->algo()->energy();
  cout << "    * CASPT2 energy:  " << setprecision(12) << setw(15) << energy_ << endl;
#endif
}


template<>
shared_ptr<GradFile> GradEval<CASPT2Grad>::compute() {
#ifdef COMPILE_SMITH
  assert(task_->target() == target_state_);
  Timer timer;

  shared_ptr<const Reference> ref = task_->ref();
  shared_ptr<FCI> fci = task_->fci();

  const int nclosed = ref->nclosed();
  const int nact = ref->nact();

  // second order density matrix
  shared_ptr<const Matrix> d1 = task_->d1();
  // first order density matrices
  shared_ptr<const Matrix> d11 = task_->d11();
  shared_ptr<const Dvec> cider = nact ? task_->cideriv() : nullptr;

  shared_ptr<const Matrix> coeff = task_->coeff();

  const int ncore = task_->ncore();
  const int nocc  = ref->nocc();
  const int nmobasis = coeff->mdim();

  // d0 including core
  auto d0ms = make_shared<Matrix>(nmobasis, nmobasis);
  if (nact)
    d0ms->add_block(1.0, nclosed, nclosed, nact, nact, task_->d10ms());
  for (int i = 0; i != nclosed; ++i)
    d0ms->element(i,i) = 2.0;

  const MatView ocoeff = coeff->slice(0, nocc);

  {
    auto dtotao = make_shared<Matrix>(*coeff * (*d0ms + *d11 + *d1) ^ *coeff);
    Dipole dipole(geom_, dtotao, "CASPT2 unrelaxed");
    dipole.compute();
  }

  // compute Yrs
  shared_ptr<const DFHalfDist> half   = ref->geom()->df()->compute_half_transform(ocoeff);
  shared_ptr<const DFHalfDist> halfj  = half->apply_J();
  shared_ptr<const DFHalfDist> halfjj = halfj->apply_J();
  shared_ptr<Matrix> yrs;
  shared_ptr<const DFFullDist> fulld1; // (gamma| ir) D(ir,js)
  tie(yrs, fulld1) = task_->compute_Y(half, halfj, halfjj);

  timer.tick_print("Yrs evaluation");

  // solve CPCASSCF
  shared_ptr<Matrix> g0 = yrs;
  shared_ptr<Dvec> g1 = nact ? cider->copy() : make_shared<Dvec>(make_shared<Determinants>(), 1);
  auto grad = make_shared<PairFile<Matrix, Dvec>>(g0, g1);

  shared_ptr<const Dvec> civector;
  if (nact) {
    civector = ref->ciwfn()->civectors();
  } else {
    auto civec = make_shared<Dvec>(make_shared<Determinants>(), 1);
    civec->data(0)->element(0,0) = 1.0;
    civector = civec;
  }
  auto cp = make_shared<CPCASSCF>(grad, civector, halfj, ref, fci, ncore, coeff);
  shared_ptr<const Matrix> zmat, xmat, smallz;
  shared_ptr<const Dvec> zvec;
  tie(zmat, zvec, xmat, smallz) = cp->solve(task_->thresh(), task_->maxziter(), task_->dcheck(), /*xms*/!!task_->dcheck());

  timer.tick_print("Z-CASSCF solution");

  // form relaxed 1RDM
  // form Zd + dZ^+
  shared_ptr<const Matrix> d0sa = nact ? ref->rdm1_mat()->resize(nmobasis, nmobasis) : d0ms;
  auto dm = make_shared<Matrix>(*zmat * *d0sa + (*d0sa ^ *zmat));

  auto dtot = make_shared<Matrix>(*d0ms + *d11 + *d1 + *dm);
  if (smallz)
    dtot->add_block(1.0, 0, 0, nocc, nocc, smallz);

  // form zdensity
  shared_ptr<const RDM<1>> zrdm1;
  shared_ptr<const RDM<2>> zrdm2;
  if (nact) {
    auto detex = make_shared<Determinants>(nact, fci->nelea(), fci->neleb(), false, /*mute=*/true);
    tie(zrdm1, zrdm2) = fci->compute_rdm12_av_from_dvec(ref->ciwfn()->civectors(), zvec, detex);

    shared_ptr<Matrix> zrdm1_mat = zrdm1->rdm1_mat(nclosed, false)->resize(nmobasis, nmobasis);
    zrdm1_mat->symmetrize();
    dtot->ax_plus_y(1.0, zrdm1_mat);
  }

  // compute relaxed dipole to check
  auto dtotao = make_shared<Matrix>(*coeff * *dtot ^ *coeff);
  {
    Dipole dipole(geom_, dtotao, "CASPT2 relaxed");
    dipole.compute();
  }

  // xmat in the AO basis
  auto xmatao = make_shared<Matrix>(*coeff * *xmat ^ *coeff);

  // two-body part
  // first make occ-occ part (copy-and-paste from src/casscf/supercigrad.cc)
  shared_ptr<const DFFullDist> qij  = halfjj->compute_second_transform(ocoeff);
  shared_ptr<DFHalfDist> qri;
  {
    shared_ptr<const Matrix> ztrans = make_shared<Matrix>(*coeff * zmat->slice(0,nocc));
    if (nact) {
      RDM<2> D(*task_->d20ms()+*zrdm2);
      RDM<1> dd(*task_->d10ms()+*zrdm1);
      // symetrize dd (zrdm1 needs symmetrization)
      for (int i = 0; i != nact; ++i)
        for (int j = 0; j != nact; ++j)
          dd(j,i) = dd(i,j) = 0.5*(dd(j,i)+dd(i,j));
      shared_ptr<DFFullDist> qijd = qij->apply_2rdm(D, dd, nclosed, nact);

      qijd->ax_plus_y(2.0, halfjj->compute_second_transform(ztrans)->apply_2rdm(*ref->rdm2_av(), *ref->rdm1_av(), nclosed, nact));
      qri = qijd->back_transform(ocoeff);

      shared_ptr<const DFFullDist> qijd2 = qij->apply_2rdm(*ref->rdm2_av(), *ref->rdm1_av(), nclosed, nact);
      qri->ax_plus_y(2.0, qijd2->back_transform(ztrans));

    } else {
      shared_ptr<DFFullDist> qijd = qij->apply_closed_2RDM();
      qijd->ax_plus_y(2.0, halfjj->compute_second_transform(ztrans)->apply_closed_2RDM());
      qri = qijd->back_transform(ocoeff);

      shared_ptr<const DFFullDist> qijd2 = qij->apply_closed_2RDM();
      qri->ax_plus_y(2.0, qijd2->back_transform(ztrans));
    }
  }

  // computing hyperfine coupling
  if (task_->do_hyperfine()) {
    shared_ptr<Matrix> dhfcc = task_->spin_density_unrelaxed();
    { // for the time being, print unrelaxed HFCC
      HyperFine hfcc(geom_, dhfcc, fci->det()->nspin(), "CASPT2 unrelaxed");
      hfcc.compute();
    }
    dhfcc->ax_plus_y(1.0, task_->spin_density_relax(zrdm1, zrdm2, zmat));
    HyperFine hfcc(geom_, dhfcc, fci->det()->nspin(), "CASPT2 relaxed");
    hfcc.compute();
  }

  // D1 part. 2.0 seems to come from the difference between smith and bagel (?)
  qri->ax_plus_y(2.0, fulld1->apply_J()->back_transform(coeff));

  // contributions from non-separable part
  shared_ptr<Matrix> qq  = qri->form_aux_2index(halfjj, 1.0);

  // separable part
  vector<shared_ptr<const Matrix>> da;
  vector<shared_ptr<const VectorB>> ca;

  auto separable_pair = [&,this](shared_ptr<const Matrix> d0occ, shared_ptr<const Matrix> d1bas) {
    shared_ptr<const Matrix> d0mo = make_shared<Matrix>(*d0occ ^ ocoeff);
    shared_ptr<const Matrix> d0ao = make_shared<Matrix>(ocoeff * *d0mo);
    shared_ptr<const Matrix> d1ao = make_shared<Matrix>(*coeff * *d1bas ^ *coeff);
    shared_ptr<const VectorB> cd0 = geom_->df()->compute_cd(d0ao);
    shared_ptr<const VectorB> cd1 = geom_->df()->compute_cd(d1ao);
    ca.push_back(cd0);
    da.push_back(d1ao);

    shared_ptr<DFHalfDist> sepd = halfjj->apply_density(d1ao);
    sepd->rotate_occ(d0occ);

    qri->ax_plus_y(-1.0, sepd);
    qri->add_direct_product(cd1, d0mo, 1.0);

    *qq += (*cd0 ^ *cd1) * 2.0;
    *qq += *halfjj->form_aux_2index(sepd, -1.0);
    return make_tuple(cd0, d1ao);
  };

  separable_pair(d0sa->get_submatrix(0,0,nocc,nocc), d1);

  if (smallz)
    separable_pair(smallz, d0sa);

  // back transform the rest
  shared_ptr<DFDist> qrs = qri->back_transform(ocoeff);
  qrs->add_direct_product(ca, da, 1.0);

  timer.tick_print("Effective densities");

  // compute gradients
  shared_ptr<GradFile> gradient = contract_gradient(dtotao, xmatao, qrs, qq);
  gradient->print();
  timer.tick_print("Gradient integral contraction");

  // set energy
  energy_ = task_->energy();
  return gradient;
#else
  return nullptr;
#endif
}


tuple<shared_ptr<Matrix>, shared_ptr<const DFFullDist>>
  CASPT2Grad::compute_Y(shared_ptr<const DFHalfDist> half, shared_ptr<const DFHalfDist> halfj, shared_ptr<const DFHalfDist> halfjj) {
#ifdef COMPILE_SMITH
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();
  const int nmobasis = coeff_->mdim();

  const MatView ocmat = coeff_->slice(0, nocc);

  shared_ptr<const DFFullDist> full = halfj->compute_second_transform(coeff_);
  shared_ptr<const DFFullDist> fullo = halfj->compute_second_transform(ocmat);

  // Y_rs = 2[Y1 + Y2 + Y3(ri) + Y4 + Y5(ri)]
  shared_ptr<Matrix> out = make_shared<Matrix>(nmobasis, nmobasis);
  const MatView ocoeff = coeff_->slice(0, nocc);

  shared_ptr<RDM<1>> rdm10;
  shared_ptr<RDM<2>> rdm20;
  {
    // 2 Y1 = h(d0 + d1 + d2) * 2
    // one-electron contributions
    const Matrix hmo(*coeff_ % *ref_->hcore() * *coeff_);
    auto d0 = make_shared<Matrix>(nmobasis, nmobasis);
    if (nact)
      d0->add_block(1.0, nclosed, nclosed, nact, nact, *d10ms_);
    for (int i = 0; i != nclosed; ++i)
      d0->element(i,i) = 2.0;
    *out += hmo * (*d1_ + *d11_ + *d0) * 2.0;
  }

  {
    // Y2 = Y2_rs = Y2_ri + Y2_ra, so making both at once
    shared_ptr<Matrix> dkl;
    if (nact) {
      dkl = ref_->rdm1_mat(); // D0SA
      dkl->sqrt();
      dkl->scale(1.0/sqrt(2.0));
    }
    Fock<1> fock(geom_, ref_->hcore()->clone(), nullptr, dkl ? ocoeff * *dkl : Matrix(ocoeff), /*grad*/false, /*rhf*/true);
    *out += *coeff_ % fock * *coeff_ * *d1_ * 2.0;
  }

  {
    // 2 Y3 = 2 Y3_ri*dm0_ji
    // coulomb
    auto dmrao = make_shared<Matrix>(*coeff_ * *d1_ ^ *coeff_);
    auto jop = geom_->df()->compute_Jop(dmrao);
    // exchange
    auto kopi = halfjj->compute_Kop_1occ(dmrao, -0.5)->transpose();

    auto tmp = make_shared<Matrix>(*coeff_ % (*jop * ocoeff + *kopi));
    if (nact) {
      *tmp *= *ref_->rdm1_mat(); // D0SA
    } else {
      *tmp *= 2.0;
    }
    out->add_block(2.0, 0, 0, nmobasis, nocc, tmp);
  }

  // fullks will be reused for gradient contraction
  shared_ptr<const DFFullDist> fullks;
  {
    // 2 Y4 =  2 K^{kl}_{rt} D^{lk}_{ts} = 2 (kr|lj) D0_(lj,ki) +  2 (kr|lt) D1_(lt,ks)
    // construct stepwise, D1 part
    fullks = contract_D1(full);
    *out += *full->form_2index(fullks, 2.0);
    // D0 part
    shared_ptr<const DFFullDist> fulld = nact ? fullo->apply_2rdm(*d20ms_, *d10ms_, nclosed, nact) : fullo->apply_closed_2RDM();
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
#else
  return tuple<shared_ptr<Matrix>, shared_ptr<const DFFullDist>>();
#endif
}

