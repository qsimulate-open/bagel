//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: casscf.cc
// Copyright (C) 2011 Toru Shiozaki
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


#include <fstream>
#include <src/scf/hf/fock.h>
#include <src/multi/casscf/casscf.h>
#include <src/multi/casscf/qvec.h>

using namespace std;
using namespace bagel;


CASSCF::CASSCF(shared_ptr<const PTree> idat, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
  : Method(idat, geom, re), hcore_(make_shared<Hcore>(geom)) {

  // drop the reference if restart is requested
  if (idata_->get<bool>("restart", false))
    ref_ = nullptr;

  if (!ref_) {
    auto scf = make_shared<RHF>(idat, geom);
    scf->compute();
    ref_ = scf->conv_to_ref();
  }

  common_init();

}


void CASSCF::common_init() {
  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("CASSCF: C1 only at the moment.");
  print_header();

  const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
  if (iactive) {
    set<int> active_indices;
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iactive) active_indices.insert(lexical_cast<int>(i->data()) - 1);
    ref_ = ref_->set_active(active_indices);
    cout << " " << endl;
    cout << "    ==== Active orbitals : ===== " << endl;
    for (auto& i : active_indices) cout << "         Orbital " << i+1 << endl;
    cout << "    ============================ " << endl << endl;
  }

  // first set coefficient
  coeff_ = ref_->coeff();
  if (geom_->nbasis() != coeff_->mdim()) {
    Overlap ovl(geom_);
    shared_ptr<const Matrix> tildex = ovl.tildex();

    Matrix c(coeff_->ndim(), tildex->mdim());
    c.copy_block(0, 0, coeff_->ndim(), coeff_->mdim(), coeff_);

    shared_ptr<const Matrix> trans = get<0>((*tildex % ovl * *coeff_).svd());
    c.copy_block(0, coeff_->mdim(), coeff_->ndim(), tildex->mdim()-coeff_->mdim(), *tildex * trans->slice(coeff_->mdim(), tildex->mdim()));
    coeff_ = make_shared<Coeff>(move(c));

#ifndef NDEBUG
    Matrix unit(coeff_->mdim(), coeff_->mdim()); unit.unit();
    assert((*coeff_ % ovl * *coeff_ - unit).rms() < 1.0e-10);
#endif
  }

  // get maxiter from the input
  max_iter_ = idata_->get<int>("maxiter", 50);
  // get maxiter from the input
  max_micro_iter_ = idata_->get<int>("maxiter_micro", 100);
  // get nstate from the input
  nstate_ = idata_->get<int>("nstate", 1);
  // get thresh (for macro iteration) from the input
  thresh_ = idata_->get<double>("thresh", 1.0e-8);
  // get thresh (for micro iteration) from the input
  thresh_micro_ = idata_->get<double>("thresh_micro", 5.0e-6);
  // whether or not to throw if the calculation does not converge
  conv_ignore_ = idata_->get<bool>("conv_ignore", false);
  // option for printing natural orbitals
  natocc_ = idata_->get<bool>("natocc", false);
  // sorting algorithm used for natural orbitals (The alternative is to sort by occupation number)
  sort_by_coeff_ = idata_->get<bool>("sort_by_coeff", true);

  // nocc from the input. If not present, full valence active space is generated.
  nact_ = idata_->get<int>("nact", 0);
  nact_ = idata_->get<int>("nact_cas", nact_);

  // nclosed from the input. If not present, full core space is generated.
  nclosed_ = idata_->get<int>("nclosed", -1);
  if (nclosed_ < -1) {
    throw runtime_error("It appears that nclosed < 0. Check nocc value.");
  } else if (nclosed_ == -1) {
    cout << "    * full core space generated for nclosed." << endl;
    nclosed_ = geom_->num_count_ncore_only() / 2;
  }
  nocc_ = nclosed_ + nact_;

  nmo_ = coeff_->mdim();
  nvirt_ = nmo_ - nocc_;
  if (nvirt_ < 0) throw runtime_error("It appears that nvirt < 0. Check the nocc value");

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl;

  const int idel = geom_->nbasis() - nmo_;
  if (idel)
    cout << "      Due to linear dependency, " << idel << (idel==1 ? " function is" : " functions are") << " omitted" << endl;


  // CASSCF methods should have FCI member. Inserting "ncore" and "norb" keyword for closed and total orbitals.
  muffle_ = make_shared<Muffle>("casscf.log");
  if (nact_) {
    auto idata = make_shared<PTree>(*idata_);
    idata->erase("active");
    fci_ = make_shared<KnowlesHandy>(idata, geom_, ref_, nclosed_, nact_, /*nstates to be read from idata*/-1, /*store*/true);
  }
  muffle_->unmute();

  do_hyperfine_ = idata_->get<bool>("hyperfine", false);

  cout <<  "  === CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

}


CASSCF::~CASSCF() {

}


void CASSCF::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "      CASSCF calculation     " << endl;
  cout << "  ---------------------------" << endl << endl;
}

void CASSCF::print_iteration(const int iter, const vector<double>& energy, const double error, const double time) const {
  muffle_->unmute();
  if (energy.size() != 1 && iter) cout << endl;

  int i = 0;
  for (auto& e : energy) {
    cout << "  " << setw(5) << iter << setw(3) << i << setw(19) << fixed << setprecision(8) << e << "   "
                 << setw(10) << scientific << setprecision(2) << (i==0 ? error : 0.0) << fixed << setw(10) << setprecision(2) << time << endl;
    ++i;
  }
  muffle_->mute();
}


shared_ptr<Matrix> CASSCF::ao_rdm1(shared_ptr<const RDM<1>> rdm1, const bool inactive_only) const {
  // first make 1RDM in MO
  const size_t nmobasis = coeff_->mdim();
  auto mo_rdm1 = make_shared<Matrix>(nmobasis, nmobasis);
  for (int i = 0; i != nclosed_; ++i) mo_rdm1->element(i,i) = 2.0;
  if (!inactive_only) {
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nact_; ++j) {
        mo_rdm1->element(nclosed_+j, nclosed_+i) = rdm1->element(j,i);
      }
    }
  }
  // transform into AO basis
  return make_shared<Matrix>(*coeff_ * *mo_rdm1 ^ *coeff_);
}


std::shared_ptr<Matrix> CASSCF::compute_active_fock(const MatView acoeff, shared_ptr<const RDM<1>> rdm1) const {
  Matrix dkl(nact_, nact_);
  copy_n(rdm1->data(), nact_*nact_, dkl.data());
  dkl.sqrt();
  dkl.scale(1.0/sqrt(2.0));
  return make_shared<Fock<1>>(geom_, hcore_->clone(), nullptr, acoeff * dkl, /*store*/false, /*rhf*/true);
}


void CASSCF::one_body_operators(shared_ptr<Matrix>& f, shared_ptr<Matrix>& fact, shared_ptr<Matrix>& factp, shared_ptr<Matrix>& gaa,
                                shared_ptr<RotFile>& d, const bool superci) const {

  shared_ptr<Matrix> finact;

  // get quantity Q_xr = 2(xs|tu)P_rs,tu (x=general)
  // note: this should be after natorb transformation.
  auto qxr = make_shared<Qvec>(coeff_->mdim(), nact_, coeff_, nclosed_, fci_, fci_->rdm2_av());

  {
    // Fock operators
    // make a matrix that contains rdm1_av
    auto rdm1mat = make_shared<Matrix>(nact_, nact_);
    copy_n(fci_->rdm1_av()->data(), rdm1mat->size(), rdm1mat->data());
    rdm1mat->sqrt();
    rdm1mat->scale(1.0/sqrt(2.0));
    auto acoeff = coeff_->slice(nclosed_, nclosed_+nact_);

    finact = make_shared<Matrix>(*coeff_ % *fci_->jop()->core_fock() * *coeff_);
    auto fact_ao = make_shared<Fock<1>>(geom_, hcore_->clone(), nullptr, acoeff * *rdm1mat, false, /*rhf*/true);
    f = make_shared<Matrix>(*finact + *coeff_% *fact_ao * *coeff_);
  }
  {
    // active-x Fock operator Dts finact_sx + Qtx
    fact = qxr->copy();// nmo_ runs first
    for (int i = 0; i != nact_; ++i)
      daxpy_(nmo_, occup_(i), finact->element_ptr(0,nclosed_+i), 1, fact->data()+i*nmo_, 1);
  }

  {
    // active Fock' operator (Fts+Fst) / (ns+nt)
    factp = make_shared<Matrix>(nact_, nact_);
    for (int i = 0; i != nact_; ++i)
      for (int j = 0; j != nact_; ++j) {
#if 1
        if (occup_(i)+occup_(j) > occup_thresh)
          factp->element(j,i) = (fact->element(j+nclosed_,i)+fact->element(i+nclosed_,j)) / (occup_(i)+occup_(j));
        else
          factp->element(j,i) = 0.0;
#else
        factp->element(j,i) = (fact->element(j+nclosed_,i)+fact->element(i+nclosed_,j)) *0.5;
#endif
      }
  }

  // G matrix (active-active) Drs,tu Factp_tu - delta_rs nr sum_v Factp_vv
  gaa = factp->clone();
  dgemv_("N", nact_*nact_, nact_*nact_, 1.0, fci_->rdm2_av()->data(), nact_*nact_, factp->data(), 1, 0.0, gaa->data(), 1);
  double p = 0.0;
  for (int i = 0; i != nact_; ++i) p += occup_(i) * factp->element(i,i);
  for (int i = 0; i != nact_; ++i) gaa->element(i,i) -= occup_(i) * p;

  // denominator
  auto denom = make_shared<RotFile>(nclosed_, nact_, nvirt_);
  fill_n(denom->data(), denom->size(), 1.0e100);

  double* target = denom->ptr_va();
  for (int i = 0; i != nact_; ++i) {
    if (occup_(i) > occup_thresh) {
      for (int j = 0; j != nvirt_; ++j, ++target)
        *target = (gaa->element(i,i) + occup_(i)*f->element(j+nocc_, j+nocc_)) / (superci ? occup_(i) : 1.0);
    } else {
      for (int j = 0; j != nvirt_; ++j, ++target)
        *target = 1.0/occup_thresh;
    }
  }

  target = denom->ptr_vc();
  for (int i = 0; i != nclosed_; ++i)
    for (int j = 0; j != nvirt_; ++j, ++target)
      *target = (f->element(j+nocc_, j+nocc_) - f->element(i, i)) / (superci ? 2.0 : 1.0);

  target = denom->ptr_ca();
  for (int i = 0; i != nact_; ++i) {
    if (2.0-occup_(i) > occup_thresh) {
      for (int j = 0; j != nclosed_; ++j, ++target)
        *target = ((f->element(nclosed_+i,nclosed_+i)*2.0-fact->element(i+nclosed_,i)) - f->element(j, j)*(2.0-occup_(i))) / (superci ? 2.0-occup_(i) : 1.0);
    } else {
      for (int j = 0; j != nclosed_; ++j, ++target)
        *target = 1.0/occup_thresh;
    }
  }
  d = denom;
}


shared_ptr<const Coeff> CASSCF::update_coeff(const shared_ptr<const Matrix> cold, shared_ptr<const Matrix> mat) const {
  auto cnew = make_shared<Coeff>(*cold);
  int nbas = cold->ndim();
  assert(nbas == geom_->nbasis());
  dgemm_("N", "N", nbas, nact_, nact_, 1.0, cold->data()+nbas*nclosed_, nbas, mat->data(), nact_,
                   0.0, cnew->data()+nbas*nclosed_, nbas);
  return cnew;
}


shared_ptr<const Coeff> CASSCF::semi_canonical_orb() const {
  auto rdm1mat = make_shared<Matrix>(nact_, nact_);
  if (nact_) {
    copy_n(fci_->rdm1_av()->data(), rdm1mat->size(), rdm1mat->data());
    rdm1mat->sqrt();
    rdm1mat->scale(1.0/sqrt(2.0));
  }
  auto ocoeff = coeff_->slice(0, nclosed_);
  auto acoeff = coeff_->slice(nclosed_, nocc_);
  auto vcoeff = coeff_->slice(nocc_, nmo_);

  VectorB eig(coeff_->mdim());
  Fock<1> fock;
  if (nact_)
    fock = Fock<1>(geom_, fci_->jop()->core_fock(), nullptr, acoeff * *rdm1mat, false, /*rhf*/true);
  else
    fock = Fock<1>(geom_, ref_->hcore(), nullptr, coeff_->slice_copy(0, nclosed_), false, /*rhf*/true);

  Matrix trans(nmo_, nmo_);
  trans.unit();
  if (nclosed_) {
    Matrix ofock = ocoeff % fock * ocoeff;
    ofock.diagonalize(eig);
    trans.copy_block(0, 0, nclosed_, nclosed_, ofock);
  }
  Matrix vfock = vcoeff % fock * vcoeff;
  vfock.diagonalize(eig);
  trans.copy_block(nocc_, nocc_, nvirt_, nvirt_, vfock);
  return make_shared<Coeff>(*coeff_ * trans);
}


shared_ptr<const Matrix> CASSCF::spin_density() const {
  Matrix den(nact_, nact_);
  shared_ptr<const RDM<1>> rdm1 = fci_->rdm1(0);
  copy_n(rdm1->data(), nact_*nact_, den.data());
  den.scale((4.0 - fci_->det()->nelea() - fci_->det()->neleb()) * 0.5);

  shared_ptr<RDM<2>> rdm2 = fci_->rdm2(0);
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      for (int k = 0; k != nact_; ++k)
        den(j,i) -= rdm2->element(j,k,k,i);

  den.scale(1.0 / (fci_->det()->nspin()*0.5 + 1.0));
  auto acoeff = coeff_->slice(nclosed_, nclosed_+nact_);
  return make_shared<Matrix>(acoeff * den ^ acoeff);
}


shared_ptr<const Reference> CASSCF::conv_to_ref() const {
  return nact_ ? make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, energy_,
                                        fci_->rdm1(), fci_->rdm2(), fci_->rdm1_av(), fci_->rdm2_av(), fci_->conv_to_ciwfn())
               : make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, energy_);
}


void CASSCF::print_natocc() const {
  assert(occup_.size() > 0);
  cout << " " << endl;
  cout << "  ========       state-averaged       ======== " << endl;
  cout << "  ======== natural occupation numbers ======== " << endl;
  for (int i=0; i!=occup_.size(); ++i)
    cout << setprecision(4) << "   Orbital " << i << " : " << occup_[i] << endl;
  cout << "  ============================================ " << endl;
}
