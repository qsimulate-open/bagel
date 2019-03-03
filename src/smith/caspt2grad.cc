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
#include <src/prop/moprint.h>


using namespace std;
using namespace bagel;

CASPT2Grad::CASPT2Grad(shared_ptr<const PTree> inp, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
  : Method(inp, geom, ref) {
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
    if (ref_->nact())
      fci_ = cas->fci();
    thresh_ = cas->thresh();
  }

  // gradient/property calculation
  do_hyperfine_ = inp->get<bool>("hyperfine", false);
  timer.tick_print("Reference calculation");

  cout << endl << "  === DF-CASPT2Grad calculation ===" << endl << endl;

#else
  throw logic_error("CASPT2 gradients require SMITH-generated code. Please compile BAGEL with --enable-smith");
#endif
}


// compute smith and set rdms and ci deriv to a member
void CASPT2Grad::compute() {
#ifdef COMPILE_SMITH

  // construct SMITH here
  shared_ptr<PTree> smithinput = idata_->get_child("smith");

  smithinput->put("_grad", true);
  smithinput->put("_hyperfine", do_hyperfine_);

  smith_ = make_shared<Smith>(smithinput, ref_->geom(), ref_);
  smith_->compute();
#endif
}


template<>
void GradEval<CASPT2Grad>::compute_dipole() const {
#ifdef COMPILE_SMITH
  const int nmobasis = ref_->coeff()->mdim();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nstate = ref_->nstate();

  vector<vector<double>> state_dipole;
  vector<vector<double>> transition_dipole;

  for (int istate = 0; istate != nstate; ++istate) {
    task_->compute_gradient(istate, istate, make_shared<NacmType>("interstate"), /*nocider=*/true);

    {
      auto d0ms = make_shared<Matrix>(nmobasis, nmobasis);
      if (nact)
        d0ms->add_block(1.0, nclosed, nclosed, nact, nact, task_->d10ms());
      for (int i = 0; i != nclosed; ++i)
        d0ms->element(i,i) = 2.0;
      {
        const string dmlabel = "CASPT2 unrelaxed dipole moment: " + to_string(istate);
        auto dtotao = make_shared<Matrix>(*(task_->smith()->algo()->coeff()) * (*d0ms + *(task_->d11()) + *(task_->d1())) ^ *(task_->smith()->algo()->coeff()));
        Dipole dipole(geom_, dtotao, dmlabel);
        auto moment = dipole.compute();
        state_dipole.push_back(moment);
      }
    }
  }

  for (int istate = 1; istate != nstate; ++istate) {
    for (int jstate = 0; jstate != istate; ++jstate) {
      task_->compute_gradient(istate, jstate, make_shared<NacmType>("interstate"), /*nocider=*/true);

      {
        auto d0ms = make_shared<Matrix>(nmobasis, nmobasis);
        if (nact)
          d0ms->add_block(1.0, nclosed, nclosed, nact, nact, task_->d10ms());
        {
          const string dmlabel = "CASPT2 unrelaxed dipole moment: " + to_string(istate) + " - " + to_string(jstate);
          auto dtotao = make_shared<Matrix>(*(task_->smith()->algo()->coeff()) * (*d0ms + *(task_->d11()) + *(task_->d1())) ^ *(task_->smith()->algo()->coeff()));
          Dipole dipole(geom_, dtotao, dmlabel);
          auto moment = dipole.compute();
          transition_dipole.push_back(moment);
        }
      }
    }
  }

  cout << "  * CASPT2 dipole moments" << endl << endl;
  for (int istate = 0; istate != nstate; ++istate) {
    cout << "    * State   " << setw(11) << istate << " : ";
    cout << "  (" << setw(12) << setprecision(6) << state_dipole[istate][0] << ", " << setw(12) << state_dipole[istate][1]
      << ", " << setw(12) << state_dipole[istate][2] << ") a.u." << endl << endl;
  }

  for (int istate = 1, counter = 0; istate != nstate; ++istate) {
    for (int jstate = 0; jstate != istate; ++jstate, ++counter) {
      cout << "    * Transition   " << setw(2) << istate << " -" << setw(2) << jstate << " : ";
        cout << "  (" << setw(12) << setprecision(6) << transition_dipole[counter][0] << ", " << setw(12) << transition_dipole[counter][1]
             << ", " << setw(12) << transition_dipole[counter][2] << ") a.u." << endl;
      const double egap = energyvec()[istate] - energyvec()[jstate];
      auto moment = transition_dipole[counter];
      const double r2 = moment[0] * moment[0] + moment[1] * moment[1] + moment[2] * moment[2];
      const double fnm = (2.0 / 3.0) * egap * r2;

      cout << "    * Oscillator strength : " << setprecision(6) << setw(10) << fnm << " a.u." << endl << endl;
    }
  }
#endif
}


void CASPT2Grad::compute_gradient(const int istate, const int jstate, shared_ptr<const NacmType> nacmtype, const bool nocider) {
#ifdef COMPILE_SMITH
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();

  target_ = istate;

  smith_->compute_gradient(istate, jstate, nacmtype, nocider);

  coeff_ = smith_->coeff();

  if (nact && !nocider)
    cideriv_ = smith_->cideriv()->copy();
  ncore_   = smith_->algo()->info()->ncore();
  wf1norm_ = smith_->wf1norm();
  msrot_   = smith_->msrot();
  nstates_ = wf1norm_.size();
  assert(msrot_->ndim() == nstates_ && msrot_->mdim() == nstates_);
  assert(nstates_ == smith_->algo()->info()->ciwfn()->nstates());

  Timer timer;

  // save correlated density matrices d(1), d(2), and ci derivatives
  auto d1tmp = make_shared<Matrix>(*smith_->dm1());
  auto d11tmp = make_shared<Matrix>(*smith_->dm11());
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
  if (smith_->dcheck()) {
    shared_ptr<const Matrix> dc = smith_->dcheck();
    assert(dc->ndim() == nact && dc->mdim() == nact);
    auto tmp = make_shared<Matrix>(nocc, nocc);
    tmp->add_block(1.0, nclosed, nclosed, nact, nact, dc);
    dcheck_ = tmp;
  }

  // spin density matrices
  if (do_hyperfine_) {
    auto sd11tmp = make_shared<Matrix>(*smith_->sdm11());
    sd11tmp->symmetrize();
    sd11tmp->scale(2.0);
    sd1_ = d1set(smith_->sdm1());
    sd11_ = d1set(sd11tmp);
  }

  if (nact) {
    d10ms_ = make_shared<RDM<1>>(nact);
    d20ms_ = make_shared<RDM<2>>(nact);
    for (int ist = 0; ist != nstates_; ++ist) {
      const double ims = msrot(ist, jstate);
      for (int jst = 0; jst != nstates_; ++jst) {
        const double jms = msrot(jst, istate);
        shared_ptr<const RDM<1>> rdm1t;
        shared_ptr<const RDM<2>> rdm2t;
        tie(rdm1t, rdm2t) = ref_->rdm12(jst, ist, /*recompute*/true);
        d10ms_->ax_plus_y(ims*jms, *rdm1t);
        d20ms_->ax_plus_y(ims*jms, *rdm2t);
      }
    }
  }

  d2_ = smith_->dm2();

  if (istate != jstate) {
    auto vd1tmp = make_shared<Matrix>(*smith_->vd1());
    vd1_ = d1set(vd1tmp);

    auto d10IJ = make_shared<Matrix>(*(ref_->rdm1_mat_tr(d10ms_)->resize(coeff_->mdim(),coeff_->mdim())));
    *vd1_ += *d10IJ;
    cout << "    * NACME Target states: " << istate << " - " << jstate << endl;

    const double energy1 = smith_->algo()->energy(istate);
    const double energy2 = smith_->algo()->energy(jstate);
    cout << "    * Energy gap is:       " << setprecision(10) << (energy1 - energy2) * au2eV__ << " eV" << endl << endl;
  } else {
    energy_ = smith_->algo()->energy(istate);
    cout << "    * CASPT2 energy:  " << setprecision(12) << setw(15) << energy_ << endl;
  }

#endif
}


template<>
vector<double> GradEval<CASPT2Grad>::energyvec() const {
#ifdef COMPILE_SMITH
  return task_->smith()->algo()->energyvec();
#else
  return vector<double>();
#endif
}


template<>
shared_ptr<GradFile> GradEval<CASPT2Grad>::compute(const string jobtitle, shared_ptr<const GradInfo> gradinfo) {
#ifdef COMPILE_SMITH
  const int istate = gradinfo->target_state();
  const int jstate = (jobtitle == "nacme") ? gradinfo->target_state2() : istate;

  task_->compute_gradient(istate, jstate, gradinfo->nacmtype(), !gradinfo->cider_eval());

  Timer timer;

  shared_ptr<const Reference> ref = task_->ref();
  auto fci = task_->fci();

  const double egap = task_->smith()->algo()->energy(jstate) - task_->smith()->algo()->energy(istate);

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
  if (jobtitle == "nacme") {
    d0ms->symmetrize();
  } else {
    for (int i = 0; i != nclosed; ++i)
      d0ms->element(i,i) = 2.0;
  }

  const MatView ocoeff = coeff->slice(0, nocc);

  if (jobtitle == "nacme") {
    const string tdmlabel = "Transition dipole moment between " + to_string(istate) + " - " + to_string(jstate);
    auto dtotao = make_shared<Matrix>(*coeff * (*d0ms + *d11 + *d1) ^ *coeff);
    Dipole dipole(geom_, dtotao, tdmlabel);
    auto moment = dipole.compute();

    const double r2 = moment[0] * moment[0] + moment[1] * moment[1] + moment[2] * moment[2];
    const double fnm = (2.0 / 3.0) * egap * r2;

    cout << "    * Oscillator strength for transition between " << istate << " - "
      << jstate << setprecision(6) << setw(10) << fnm << " a.u." << endl << endl;
  } else {
    auto dtotao = make_shared<Matrix>(*coeff * (*d0ms + *d11 + *d1) ^ *coeff);
    Dipole dipole(geom_, dtotao, "CASPT2 unrelaxed");
    dipole.compute();
  }

  if (!gradinfo->cider_eval()) {
    return make_shared<GradFile>(geom_->natom());
  }

  // compute Yrs
  shared_ptr<const DFHalfDist> half   = ref->geom()->df()->compute_half_transform(ocoeff);
  shared_ptr<const DFHalfDist> halfj  = half->apply_J();
  shared_ptr<const DFHalfDist> halfjj = halfj->apply_J();
  shared_ptr<Matrix> yrs;
  shared_ptr<const DFFullDist> fulld1; // (gamma| ir) D(ir,js)
  tie(yrs, fulld1) = task_->compute_Y(half, halfj, halfjj, /*nacme=*/(jobtitle=="nacme"));

  timer.tick_print("Yrs evaluation");

  // solve CPCASSCF
  shared_ptr<Matrix> g0 = yrs;
  shared_ptr<Dvec> g1 = nact ? cider->copy() : make_shared<Dvec>(make_shared<Determinants>(), 1);

  if (jobtitle == "nacme" && (gradinfo->nacmtype()->is_full() || gradinfo->nacmtype()->is_etf()))
    task_->augment_Y(d0ms, g0, g1, halfj, istate, jstate, egap);

  timer.tick_print("Yrs non-Lagrangian terms");

  auto grad = make_shared<PairFile<Matrix, Dvec>>(g0, g1);

  shared_ptr<const Dvec> civector;
  if (nact) {
    civector = ref->ciwfn()->civectors();
  } else {
    auto civec = make_shared<Dvec>(make_shared<Determinants>(), 1);
    civec->data(0)->element(0,0) = 1.0;
    civector = civec;
  }
  auto cp = make_shared<CPCASSCF>(grad, civector, halfj, ref, fci, ncore, task_->smith()->algo()->info()->shift_imag(), coeff);
  shared_ptr<const Matrix> zmat, xmat, smallz;
  shared_ptr<const Dvec> zvec;
  tie(zmat, zvec, xmat, smallz) = cp->solve(task_->thresh(), gradinfo->maxziter(), task_->dcheck(), /*xms*/!!task_->dcheck());

  timer.tick_print("Z-CASSCF solution");

  // form relaxed 1RDM
  // form Zd + dZ^+
  shared_ptr<const Matrix> d0sa = nact ? ref->rdm1_mat()->resize(nmobasis, nmobasis) : d0ms;
  auto dm = make_shared<Matrix>(*zmat * *d0sa + (*d0sa ^ *zmat));
  auto dtot = make_shared<Matrix>(*d0ms + *d11 + *d1 + *dm);
  if (smallz)
    *dtot += *smallz;

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
    dipole_ = dipole.compute();
  }

  // print relaxed density if requested
  if (gradinfo->density_print()) {
    auto density_print = make_shared<MOPrint>(gradinfo->moprint_info(), geom_, ref_, /*is_density=*/true, make_shared<const ZMatrix>(*dtotao, 1.0));
    density_print->compute();
  }

  // xmat in the AO basis
  auto xmatao = make_shared<Matrix>(*coeff * *xmat ^ *coeff);
  shared_ptr<Matrix> qxmatao;
  if (jobtitle == "nacme") {
    auto qxmat = task_->vd1()->resize(nmobasis, nmobasis);

    if (gradinfo->nacmtype()->is_full())
      qxmat->scale(egap);
    else
      qxmat->zero();

    qxmatao = make_shared<Matrix>(*coeff * *qxmat ^ *coeff);
  }

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
      shared_ptr<DFFullDist> qijd;
      if (jobtitle == "nacme")
        qijd = qij->apply_2rdm_tran(D, dd, nclosed, nact);
      else
        qijd = qij->apply_2rdm(D, dd, nclosed, nact);

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

  // separable part: see Gyorffy appendix
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
    sepd = sepd->transform_occ(d0occ);

    qri->ax_plus_y(-1.0, sepd);
    qri->add_direct_product(cd1, d0mo, 1.0);

    *qq += (*cd0 ^ *cd1) * 2.0;
    *qq += *halfjj->form_aux_2index(sepd, -1.0);
    return make_tuple(cd0, d1ao);
  };

  separable_pair(d0sa->get_submatrix(0,0,nocc,nocc), d1);

  if (smallz)
    separable_pair(d0sa->get_submatrix(0,0,nocc,nocc), smallz);

  // back transform the rest
  shared_ptr<DFDist> qrs = qri->back_transform(ocoeff);
  qrs->add_direct_product(ca, da, 1.0);

  timer.tick_print("Effective densities");

  // compute gradients
  shared_ptr<GradFile> gradient = contract_gradient(dtotao, xmatao, qrs, qq, qxmatao);

  if ((jobtitle == "nacme") && !(gradinfo->nacmtype()->is_noweight()))
    gradient->scale(1.0/egap);
  gradient->print();
  timer.tick_print("Gradient integral contraction");

  if (jobtitle == "force")
    energy_ = task_->energy();
  else
    energy_ = 0.0;

  return gradient;
#else
  throw logic_error("CASPT2 gradients require SMITH-generated code. Please compile BAGEL with --enable-smith");
  return nullptr;
#endif
}
