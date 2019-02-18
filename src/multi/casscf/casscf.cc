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
  : Method(idat, geom, re), hcore_(make_shared<Hcore>(geom, geom->hcoreinfo())) {

  // check if RDMs are supplied externally
  external_rdm_ = idata_->get<string>("external_rdm", "");
  if (!external_rdm_.empty() && external_rdm_ != "noref") {
    IArchive ar("ref");
    ar >> ref_;
  }

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

  thresh_overlap_ = idata_->get<double>("thresh_overlap", 1.0e-8);

  // first set coefficient
  coeff_ = ref_->coeff();
  if (geom_->nbasis() != coeff_->mdim()) {
    Overlap ovl(geom_);
    shared_ptr<const Matrix> tildex = ovl.tildex(thresh_overlap_);

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
  energy_.resize(nstate_);
  // get thresh (for macro iteration) from the input
  thresh_ = idata_->get<double>("thresh", 1.0e-8);
  // get thresh (for micro iteration) from the input
  thresh_micro_ = idata_->get<double>("thresh_micro", 5.0e-6);

  // whether or not to throw if the calculation does not converge
  conv_ignore_ = idata_->get<bool>("conv_ignore", false);

  // to save binary archives with each iteration
  restart_cas_ = idata_->get<bool>("restart_cas", false);

  // option for printing natural orbitals
  natocc_ = idata_->get<bool>("natocc", false);

  // option for saving canonical orbitals
  canonical_ = idata_->get<bool>("canonical", false);

  // nact from the input.
  nact_ = idata_->get<int>("nact", 0);
  // FCI algorithm
  auto fci_algo = idata_->get<string>("fci_algorithm", ((nact_ > 9) && (mpi__->size() >= 8)) ? "parallel" : "knowles");
  fci_algorithm_ = make_shared<FCI_algorithms>(fci_algo);

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
    if (fci_algorithm_->is_knowles()) {
      cout << "    * Using serial Knowles-Handy algorithm in FCI." << endl;
      fci_ = make_shared<KnowlesHandy>(idata, geom_, ref_, nclosed_, nact_, /*nstates to be read from idata*/-1, /*store*/true);
    } else if (fci_algorithm_->is_harrison()) {
      cout << "    * Using serial Harrison-Zarrabian algorithm in FCI." << endl;
      fci_ = make_shared<HarrisonZarrabian>(idata, geom_, ref_, nclosed_, nact_, /*nstates to be read from idata*/-1, /*store*/true);
#ifdef HAVE_MPI_H
    } else if (fci_algorithm_->is_dist()) {
      cout << "    * Using parallel algorithm in FCI." << endl;
      fci_ = make_shared<DistFCI>(idata, geom_, ref_, nclosed_, nact_, /*nstates to be read from idata*/-1, /*store*/true);
#endif
    }
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


tuple<shared_ptr<const Coeff>,VectorB,VectorB> CASSCF::semi_canonical_orb() const {
  auto rdm1mat = make_shared<Matrix>(nact_, nact_);
  shared_ptr<const Matrix> rdm1copy; 
  if (nact_) {
    copy_n(fci_->rdm1_av()->data(), rdm1mat->size(), rdm1mat->data());
    rdm1copy = rdm1mat->copy();
    rdm1mat->sqrt();
    rdm1mat->scale(1.0/sqrt(2.0));
  }
  auto ocoeff = coeff_->slice(0, nclosed_);
  auto acoeff = coeff_->slice(nclosed_, nocc_);
  auto vcoeff = coeff_->slice(nocc_, nmo_);

  Fock<1> fock;
  if (nact_)
    fock = Fock<1>(geom_, fci_->jop()->core_fock(), nullptr, acoeff * *rdm1mat, false, /*rhf*/true);
  else
    fock = Fock<1>(geom_, ref_->hcore(), nullptr, coeff_->slice_copy(0, nclosed_), false, /*rhf*/true);

  VectorB eig(coeff_->mdim());
  VectorB occup(coeff_->mdim());
  // calculate the transformation matrix to (semi-)canonical orbitals
  Matrix trans(nmo_, nmo_);
  trans.unit();
  if (nclosed_) {
    Matrix ofock = ocoeff % fock * ocoeff;
    VectorB tmp(nclosed_);
    ofock.diagonalize(tmp);
    trans.copy_block(0, 0, nclosed_, nclosed_, ofock);
    copy_n(tmp.data(), nclosed_, eig.data());
    fill_n(occup.data(), nclosed_, 2.0);
  }
  if (nact_ && canonical_) {
    Matrix afock = acoeff % fock * acoeff;
    VectorB tmp(nact_);
    afock.diagonalize(tmp);
    trans.copy_block(nclosed_, nclosed_, nact_, nact_, afock);
    copy_n(tmp.data(), nact_, eig.data()+nclosed_); 
  }
  {
    Matrix vfock = vcoeff % fock * vcoeff;
    VectorB tmp(nvirt_);
    vfock.diagonalize(tmp);
    trans.copy_block(nocc_, nocc_, nvirt_, nvirt_, vfock);
    copy_n(tmp.data(), nvirt_, eig.data()+nclosed_+nact_);
  }
  // finally calculate the occupation numbers for active orbitals (diagonal elements of transformed 1RDM)  
  {
    shared_ptr<const Matrix> atrans = trans.get_submatrix(nclosed_, nclosed_, nact_, nact_);
    const Matrix transrdm = *atrans * *rdm1copy ^ *atrans; 
    for (int i = 0; i != nact_; ++i)
      occup[i+nclosed_] = transrdm(i, i);
  }

  auto coeffout = make_shared<Coeff>(*coeff_ * trans);
  return make_tuple(coeffout, move(eig), move(occup));
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
  const bool noci = !nact_ || external_rdm_ == "noref";
  auto ref = noci ? make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, energy_)
                  : make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, energy_,
                                           fci_->rdm1(), fci_->rdm2(), fci_->rdm1_av(), fci_->rdm2_av(), fci_->conv_to_ciwfn());
  ref->set_eig(eig_);
  ref->set_occup(occup_);
  return ref;
}
