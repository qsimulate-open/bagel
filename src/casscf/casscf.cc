//
// BAGEL - Parallel electron correlation program.
// Filename: casscf.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <src/casscf/casscf.h>
#include <src/casscf/qvec.h>

using namespace std;
using namespace bagel;

template<typename T>
static string tostring(const T i) {
  stringstream ss;
  ss << i;
  return ss.str();
};

CASSCF::CASSCF(multimap<string, string> idat, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
  : idata_(idat), geom_(geom), hcore_(new Hcore(geom)) {

  std::shared_ptr<SCF<1>> scf_;
  scf_ = std::shared_ptr<SCF<1>>(new SCF<1>(idat, geom, re));
  scf_->compute();
  ref_ = scf_->conv_to_ref();

  common_init();

}


void CASSCF::common_init() {
  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("CASSCF: C1 only at the moment.");
  print_header();

  // first set coefficient
  coeff_ = ref_->coeff();

  // get maxiter from the input
  max_iter_ = read_input<int>(idata_, "maxiter", 100);
  // get maxiter from the input
  max_micro_iter_ = read_input<int>(idata_, "maxiter_micro", 100);
  // get nstate from the input
  nstate_ = read_input<int>(idata_, "nstate", 1);
  // get istate from the input (for geometry optimization)
  istate_ = read_input<int>(idata_, "istate", 0);
  // get thresh (for macro iteration) from the input
  thresh_ = read_input<double>(idata_, "thresh", 1.0e-10);
  // get thresh (for micro iteration) from the input
  thresh_micro_ = read_input<double>(idata_, "thresh_micro", thresh_);

  // nocc from the input. If not present, full valence active space is generated.
  nact_ = read_input<int>(idata_, "nact", 0);
  nact_ = read_input<int>(idata_, "nact_cas", nocc_);

  // nclosed from the input. If not present, full core space is generated.
  nclosed_ = read_input<int>(idata_, "nclosed", -1);
  if (nclosed_ < -1) {
    throw runtime_error("It appears that nclosed < 0. Check nocc value.");
  } else if (nclosed_ == -1) {
    cout << "    * full core space generated for nclosed." << endl;
    nclosed_ = geom_->num_count_ncore_only() / 2;
  }
  nocc_ = nclosed_ + nact_;

#if 0
  nbasis_ = geom_->nbasis();
#else
  nbasis_ = coeff_->mdim();
#endif
  nvirt_ = nbasis_ - nocc_;
  if (nvirt_ < 0) throw runtime_error("It appears that nvirt < 0. Check the nocc value");

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl;

  const int idel = geom_->nbasis() - nbasis_;
  if (idel)
    cout << "      Due to linear dependency, " << idel << (idel==1 ? " function is" : " functions are") << " omitted" << endl;


  // CASSCF methods should have FCI member. Inserting "ncore" and "norb" keyword for closed and total orbitals.
  mute_stdcout();
  fci_ = shared_ptr<FCI>(new KnowlesHandy(idata_, ref_, nclosed_, nact_)); // nstate does not need to be specified as it is in idata_...
  resume_stdcout();


  schwarz_ = geom_->schwarz();

  cout <<  "  === CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

};


CASSCF::~CASSCF() {

};


void CASSCF::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "      CASSCF calculation     " << endl;
  cout << "  ---------------------------" << endl << endl;
}

void CASSCF::print_iteration(int iter, int miter, int tcount, const vector<double> energy, const double error, const double time) const {
  if (energy.size() != 1 && iter) cout << endl;
  int i = 0;
  for (auto eiter = energy.begin(); eiter != energy.end(); ++eiter, ++i)
  cout << "  " << setw(5) << iter << setw(3) << i << setw(4) << miter << setw(4) << tcount
               << setw(20) << fixed << setprecision(12) << *eiter << "   "
               << setw(10) << scientific << setprecision(2) << (i==0 ? error : 0.0) << fixed << setw(10) << setprecision(2)
               << time << endl;
}

static streambuf* backup_stream_;
static ofstream* ofs_;

void CASSCF::mute_stdcout() {
#if 1
  ofstream* ofs(new ofstream("casscf.log",(backup_stream_ ? ios::app : ios::trunc)));
  ofs_ = ofs;
  backup_stream_ = cout.rdbuf(ofs->rdbuf());
#endif
}


void CASSCF::resume_stdcout() {
#if 1
  cout.rdbuf(backup_stream_);
  delete ofs_;
#endif
}


shared_ptr<Matrix> CASSCF::ao_rdm1(shared_ptr<RDM<1>> rdm1, const bool inactive_only) const {
  // first make 1RDM in MO
  shared_ptr<Matrix> mo_rdm1(new Matrix(geom_->nbasis(), geom_->nbasis()));
  for (int i = 0; i != nclosed_; ++i) mo_rdm1->element(i,i) = 2.0;
  if (!inactive_only) {
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nact_; ++j) {
        mo_rdm1->element(nclosed_+j, nclosed_+i) = rdm1->element(j,i);
      }
    }
  }
  // transform into AO basis
  return shared_ptr<Matrix>(new Matrix(*coeff_ * *mo_rdm1 ^ *coeff_));
}



void CASSCF::one_body_operators(shared_ptr<Matrix>& f, shared_ptr<Matrix>& fact, shared_ptr<Matrix>& factp, shared_ptr<Matrix>& gaa,
                                shared_ptr<RotFile>& d, const bool superci) const {

  shared_ptr<Matrix> finact;

  // get quantity Q_xr = 2(xs|tu)P_rs,tu (x=general)
  // note: this should be after natorb transformation.
  shared_ptr<Qvec> qxr(new Qvec(geom_->nbasis(), nact_, geom_->df(), coeff_, nclosed_, fci_, fci_->rdm2_av()));

  {
    // Fock operators
    if (nclosed_) {
      shared_ptr<Matrix> deninact = ao_rdm1(fci_->rdm1_av(), true); // true means inactive_only
      finact = shared_ptr<Matrix>(new Matrix(*coeff_ % *fci_->jop()->core_fock() * *coeff_));

      shared_ptr<Matrix> denall = ao_rdm1(fci_->rdm1_av());
      shared_ptr<Matrix> denact(new Matrix(*denall-*deninact));
      shared_ptr<Fock<1>> fact_ao(new Fock<1>(geom_, hcore_, denact, schwarz_));
      f = shared_ptr<Matrix>(new Matrix(*finact + *coeff_%(*fact_ao-*hcore_)**coeff_));
    } else {
      shared_ptr<Matrix> denall = ao_rdm1(fci_->rdm1_av());
      shared_ptr<Fock<1>> f_ao(new Fock<1>(geom_, hcore_, denall, schwarz_));
      f = shared_ptr<Matrix>(new Matrix(*coeff_ % *f_ao * *coeff_));

      finact = shared_ptr<Matrix>(new Matrix(*coeff_ % *hcore_ * *coeff_));
    }
  }
  {
    // active-x Fock operator Dts finact_sx + Qtx
    fact = shared_ptr<Matrix>(new Matrix(*qxr));// nbasis_ runs first
    for (int i = 0; i != nact_; ++i)
      daxpy_(nbasis_, occup_[i], finact->element_ptr(0,nclosed_+i), 1, fact->data()+i*nbasis_, 1);
  }

  {
    // active Fock' operator (Fts+Fst) / (ns+nt)
    factp = shared_ptr<Matrix>(new Matrix(nact_, nact_));
    for (int i = 0; i != nact_; ++i)
      for (int j = 0; j != nact_; ++j) {
#if 1
        if (occup_[i]+occup_[j] > occup_thresh)
          factp->element(j,i) = (fact->element(j+nclosed_,i)+fact->element(i+nclosed_,j)) / (occup_[i]+occup_[j]);
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
  for (int i = 0; i != nact_; ++i) p += occup_[i] * factp->element(i,i);
  for (int i = 0; i != nact_; ++i) gaa->element(i,i) -= occup_[i] * p;

  // denominator
  shared_ptr<RotFile> denom(new RotFile(nclosed_, nact_, nvirt_));
  fill(denom->data(), denom->data()+denom->size(), 1.0e100);

  double* target = denom->ptr_va();
  for (int i = 0; i != nact_; ++i) {
    if (occup_[i] > occup_thresh) {
      for (int j = 0; j != nvirt_; ++j, ++target)
        *target = (gaa->element(i,i) + occup_[i]*f->element(j+nocc_, j+nocc_)) / (superci ? occup_[i] : 1.0);
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
    if (2.0-occup_[i] > occup_thresh) {
      for (int j = 0; j != nclosed_; ++j, ++target)
        *target = ((f->element(nclosed_+i,nclosed_+i)*2.0-fact->element(i+nclosed_,i)) - f->element(j, j)*(2.0-occup_[i])) / (superci ? 2.0-occup_[i] : 1.0);
    } else {
      for (int j = 0; j != nclosed_; ++j, ++target)
        *target = 1.0/occup_thresh;
    }
  }
  d = denom;
}


shared_ptr<const Coeff> CASSCF::update_coeff(const shared_ptr<const Coeff> cold, shared_ptr<Matrix> mat) const {
  shared_ptr<const Matrix> cnew(new const Matrix(*dynamic_cast<const Matrix*>(cold.get())));
  int nbas = geom_->nbasis();
  dgemm_("N", "N", nbas, nact_, nact_, 1.0, cold->data()+nbas*nclosed_, nbas, mat->data(), nact_,
                   0.0, cnew->data()+nbas*nclosed_, nbas);
  return shared_ptr<const Coeff>(new Coeff(*cnew));
}



shared_ptr<Matrix> CASSCF::form_natural_orbs() {
    // here make a natural orbitals and update the coefficients
    // this effectively updates 1,2RDM and integrals
    const pair<shared_ptr<Matrix>, vector<double>> natorb = fci_->natorb_convert();
    // new coefficients
    shared_ptr<const Coeff> new_coeff = update_coeff(coeff_, natorb.first);
    coeff_ = new_coeff;
    // occupation number of the natural orbitals
    occup_ = natorb.second;
    return natorb.first;
}


shared_ptr<const Reference> CASSCF::conv_to_ref() const {
  shared_ptr<Reference> out(new Reference(geom_, coeff_, nclosed_, nact_, nvirt_, energy(),
                                          fci_->rdm1(), fci_->rdm2(), fci_->rdm1_av(), fci_->rdm2_av()));

  // TODO
  // compute one-boedy operators
  shared_ptr<Matrix> f;
  shared_ptr<Matrix> fact, factp, gaa;
  shared_ptr<RotFile>  denom;
  one_body_operators(f, fact, factp, gaa, denom);

  *f *= 2.0;

  for (int i = 0; i != nbasis_; ++i) {
    for (int j = 0; j != nbasis_; ++j) {
      if (i < nocc_ && j < nocc_) continue;
      f->element(j,i) = 0.0;
    }
  }
  for (int j = 0; j != nact_; ++j) {
    for (int i = 0; i != nocc_; ++i) {
      f->element(i,j+nclosed_) = fact->element(i,j);
    }
  }

  shared_ptr<Matrix> erdm(new Matrix(*coeff_ * *f ^ *coeff_));

  out->set_erdm1(erdm);
  out->set_nstate(nstate_);
  return out;
}
