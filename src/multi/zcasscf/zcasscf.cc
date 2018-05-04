//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcasscf.cc
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

#include <fstream>
#include <src/scf/dhf/dirac.h>
#include <src/scf/dhf/dfock.h>
#include <src/util/math/quatmatrix.h>
#include <src/ci/zfci/relfci.h>
#include <src/multi/zcasscf/zcasscf.h>
#include <src/mat1e/giao/relhcore_london.h>
#include <src/mat1e/giao/reloverlap_london.h>
#include <src/util/atommap.h>

using namespace std;
using namespace bagel;

ZCASSCF::ZCASSCF(shared_ptr<const PTree> idat, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
  : Method(idat, geom, ref) {
  // check if RDMs are supplied externally
  external_rdm_ = idata_->get<string>("external_rdm", "");
  if (!external_rdm_.empty() && external_rdm_ != "noref") {
    IArchive ar("relref");
    shared_ptr<RelReference> r;
    ar >> r;
    ref_ = ref = r;
  }

  if (!dynamic_pointer_cast<const RelReference>(ref)) {
    if (ref != nullptr && ref->coeff()->ndim() == geom->nbasis()) {
      nr_coeff_ = ref->coeff();
      if (geom->magnetism()) // TODO Implement nr_coeff_ for GIAO basis sets
        throw runtime_error("So far only relativistic input coefficients can be used for ZCASSCF with magnetic field");
    }
  }

  // relref needed for many things below ; TODO eliminate dependence on ref_ being a relref
  if (!dynamic_pointer_cast<const RelReference>(ref)) {
    auto idata_tmp = make_shared<PTree>(*idata_);
    const int ctmp = idata_->get<int>("charge", 0);
    const int nele = geom->nele();
    if ((nele - ctmp) % 2) {
      idata_tmp->erase("charge");
      idata_tmp->put<int>("charge", ctmp - 1);
    }

    // With a CASSCF initial guess, the DHF wavefunction will not be used...
    if (nr_coeff_)
      idata_tmp->put<double>("thresh", 1.0e-4);

    auto scf = make_shared<Dirac>(idata_tmp, geom_, ref);
    scf->compute();
    ref_ = scf->conv_to_ref();
  }

}


void ZCASSCF::init() {
  print_header();

  auto relref = dynamic_pointer_cast<const RelReference>(ref_);

  gaunt_ = idata_->get<bool>("gaunt", relref->gaunt());
  breit_ = idata_->get<bool>("breit", relref->breit());

  if (!geom_->dfs() || (gaunt_ != relref->gaunt()))
    geom_ = geom_->relativistic(gaunt_);

  nneg_ = geom_->nbasis()*2;

  // nact from the input.
  nact_ = idata_->get<int>("nact", 0);
  if (!nact_) energy_.resize(1);
  // option for printing natural orbital occupation numbers
  natocc_ = idata_->get<bool>("natocc", false);

  // option for saving canonical orbitals
  canonical_ = idata_->get<bool>("canonical", false);

  // nclosed from the input. If not present, full core space is generated.
  nclosed_ = idata_->get<int>("nclosed", -1);
  if (nclosed_ < -1) {
    throw runtime_error("It appears that nclosed < 0. Check nocc value.");
  } else if (nclosed_ == -1) {
    cout << "    * full core space generated for nclosed." << endl;
    nclosed_ = geom_->num_count_ncore_only() / 2;
  }
  nocc_ = nclosed_ + nact_;

  nbasis_ = geom_->nbasis()*2;
  nvirt_ = nbasis_ - nocc_;
  if (nvirt_ < 0) throw runtime_error("It appears that nvirt < 0. Check the nocc value");
  nvirtnr_ = nvirt_ - nneg_/2;

  charge_ = idata_->get<int>("charge", 0);
  if (nclosed_*2 > geom_->nele() - charge_)
    throw runtime_error("too many closed orbitals in the input");
  thresh_overlap_ = idata_->get<double>("thresh_overlap", 1.0e-8);

  // get maxiter from the input
  max_iter_ = idata_->get<int>("maxiter", 100);
  // get maxiter from the input
  max_micro_iter_ = idata_->get<int>("maxiter_micro", 20);

  // whether or not to throw if the calculation does not converge
  conv_ignore_ = idata_->get<bool>("conv_ignore", false);

  // to save binary archives with each iteration
  restart_cas_ = idata_->get<bool>("restart_cas", false);

  batchsize_ = idata_->get<int>("batchsize", 250);
  if (batchsize_ < 2*nclosed_) {
    const int nbatch = (2*nclosed_ - 1) / batchsize_ + 1;
    const string sizerange = ((2*nclosed_) % nbatch == 0) ? to_string((2*nclosed_) / nbatch) : to_string((2*nclosed_) / nbatch) + " to " + to_string((2*nclosed_) / nbatch + 1);
    cout << "   ***  Maximum batchsize set to " << batchsize_ << ".  The core Fock matrix will be evaluated in " << nbatch << " batches of " << sizerange << " spin-orbitals." << endl;
  }

  // get thresh (for macro iteration) from the input
  thresh_ = idata_->get<double>("thresh", 1.0e-8);
  // get thresh (for micro iteration) from the input
  thresh_micro_ = idata_->get<double>("thresh_micro", thresh_*0.5);

  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl;

  cout << "    * gaunt    : " << (gaunt_ ? "true" : "false") << endl;
  cout << "    * breit    : " << (breit_ ? "true" : "false") << endl;
  cout << "    * active space: " << geom_->nele() - charge_ - nclosed_*2 << " electrons in " << nact_ << " orbitals" << endl;

  // initialize hcore and overlap
  init_mat1e();

  // initialize coefficient
  init_coeff();
  // select active orbitals if "active" is set
  select_active();

  const int idel = geom_->nbasis()*2 - nbasis_;
  if (idel)
    cout << "      Due to linear dependency, " << idel << (idel==1 ? " function is" : " functions are") << " omitted" << endl;

  muffle_ = make_shared<Muffle>("casscf.log");
  // CASSCF methods should have FCI member. Inserting "ncore" and "norb" keyword for closed and active orbitals.
  if (nact_)
    fci_ = make_shared<RelFCI>(idata_, geom_, ref_, nclosed_, nact_, coeff_, /*store*/true);
  nstate_ = nact_ ? fci_->nstate() : 1;
  energy_.resize(nstate_);
  muffle_->unmute();
  cout << "    * nstate   : " << setw(6) << nstate_ << endl << endl;

  cout <<  "  === Dirac CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

}


void ZCASSCF::select_active() {
  // specify active orbitals and move into the active space
  shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
  if (iactive) {
    shared_ptr<const ZCoeff_Striped> scoeff = coeff_->striped_format();
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    set<int> active_indices;
    for (auto& i : *iactive)
      active_indices.insert(lexical_cast<int>(i->data()) - 1);
    scoeff = scoeff->set_active(active_indices, geom_->nele()-charge_);
    coeff_ = scoeff->block_format();
  }
}


void ZCASSCF::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "      CASSCF calculation     " << endl;
  cout << "  ---------------------------" << endl << endl;
}


void ZCASSCF::print_iteration(const int iter, const vector<double>& energy, const double error, const double time) const {
  muffle_->unmute();
  if (energy.size() != 1 && iter) cout << endl;
  int i = 0;
  for (auto& e : energy) {
    cout << "     " << setw(5) << iter << setw(4) << i << " " << setw(19) << fixed << setprecision(8) << e << "   "
                    << setw(10) << scientific << setprecision(2) << (i==0 ? error : 0.0) << fixed << setw(10) << setprecision(2) << time << endl;
    ++i;
  }
  muffle_->mute();
}


shared_ptr<const Reference> ZCASSCF::conv_to_ref_(const bool kramers) const {
  const bool noci = !nact_ || external_rdm_ == "noref";
  return noci ? make_shared<RelReference>(geom_, coeff_->striped_format(), energy_, nneg_, nclosed_, nact_, nvirtnr_, gaunt_, breit_, kramers)
              : make_shared<RelReference>(geom_, coeff_->striped_format(), energy_, nneg_, nclosed_, nact_, nvirtnr_, gaunt_, breit_, kramers,
                                          fci_->rdm1_av(), fci_->rdm2_av(), fci_->conv_to_ciwfn());
}


shared_ptr<ZMatrix> ZCASSCF::compute_active_fock(const ZMatView coeff, shared_ptr<const ZMatrix> rdm1, const bool coulomb_only) const {
  // calculate S^1/2 of rdm1
  ZMatrix s(*rdm1);
  s.sqrt();
  const bool do_gaunt = !coulomb_only && gaunt_;
  const bool do_breit = !coulomb_only && breit_;
  auto zero = make_shared<ZMatrix>(coeff.ndim(), coeff.ndim());
  return make_shared<DFock>(geom_, zero, coeff * *s.get_conjg(), do_gaunt, do_breit, /*store half*/false, /*robust*/do_breit);
}


shared_ptr<ZCoeff_Kramers> ZCASSCF::nonrel_to_relcoeff(shared_ptr<const Matrix> nr_coeff) const {
  // constructs a relativistic coefficient for electronic components from a non-rel coefficient
  const int n = nr_coeff->ndim();
  const int m = nr_coeff->mdim();
  assert(nvirt_ - nneg_/2 == nvirtnr_);

  // compute T^(-1/2)
  shared_ptr<ZMatrix> t12 = overlap_->get_submatrix(n*2, n*2, n, n);
  t12 = t12->tildex(thresh_overlap_);

  // compute S^(1/2)
  shared_ptr<ZMatrix> shalf = overlap_->get_submatrix(0, 0, n, n);
  shalf = shalf->tildex(thresh_overlap_);

  if (t12->mdim() != shalf->mdim())
    throw runtime_error("Different linear dependency for the overlap and kinetic matrices in conversion to relativistic coefficients.");

  // compute positronic orbital coefficients
  auto tcoeff = make_shared<ZMatrix>(n, m);
  tcoeff->add_real_block(1.0, 0, 0, n, m, *nr_coeff);
  *tcoeff = *t12 * (*shalf % *tcoeff);

  // build output coefficient matrix
  auto out = make_shared<ZCoeff_Kramers>(4*n, nr_coeff->localized(), nclosed_, nact_, nvirtnr_, nneg_);
  assert(out->mdim() == 4*tcoeff->mdim());
  out->copy_real_block(1.0, 0, 0, n, m, *nr_coeff);
  out->copy_real_block(1.0, n, 2*m, n, m, *nr_coeff);
  out->copy_block(2*n, m, n, m, *tcoeff);
  out->copy_block(3*n, 3*m, n, m, *tcoeff);
  return out;
}


// Eliminates the positronic entries for the given rot file
void ZCASSCF::zero_positronic_elements(shared_ptr<ZRotFile> rot) {
  int nr_nvirt = nvirt_ - nneg_/2;
  for (int i = 0; i != nclosed_*2; ++i) {
    for (int j = 0; j != nneg_/2; ++j) {
      rot->ele_vc(j + nr_nvirt, i) = 0.0;
      rot->ele_vc(j + nr_nvirt + nvirt_, i) = 0.0;
    }
  }
  for (int i = 0; i != nact_*2; ++i) {
    for (int j = 0; j != nneg_/2; ++j) {
      rot->ele_va(j + nr_nvirt, i) = 0.0;
      rot->ele_va(j + nr_nvirt + nvirt_, i) = 0.0;
    }
  }
}


shared_ptr<const ZCoeff_Block> ZCASSCF::semi_canonical_orb(const bool kramers) const {
  const ZMatrix fock = *fci_->jop()->core_fock()
                     + *compute_active_fock(coeff_->slice(nclosed_*2, nocc_*2), fci_->rdm1_av());
  auto ocoeff = coeff_->slice(0, nclosed_*2);
  auto acoeff = coeff_->slice(nclosed_*2, nocc_*2);
  auto vcoeff = coeff_->slice(nocc_*2, (nocc_+nvirt_)*2);
  assert((nocc_+nvirt_)*2 == coeff_->mdim());

  ZMatrix trans(coeff_->mdim(), coeff_->mdim());
  trans.unit();
  VectorB eig(coeff_->mdim());

  if (kramers) {
    if (nclosed_) {
      QuatMatrix ofock(ocoeff % fock * ocoeff);
      ofock.diagonalize(eig);
      trans.copy_block(0, 0, nclosed_*2, nclosed_*2, ofock);
    }
    if (nact_ && canonical_) {
      QuatMatrix afock = acoeff % fock * acoeff;
      afock.diagonalize(eig);
      trans.copy_block(nclosed_*2, nclosed_*2, nact_*2, nact_*2, afock);
    }
    QuatMatrix vfock = vcoeff % fock * vcoeff;
    vfock.diagonalize(eig);
    trans.copy_block(nocc_*2, nocc_*2, nvirt_*2, nvirtnr_, vfock.slice(nvirt_-nvirtnr_, nvirt_));
    trans.copy_block(nocc_*2, nocc_*2+nvirtnr_, nvirt_*2, nvirt_-nvirtnr_, vfock.slice(0, nvirt_-nvirtnr_));
    trans.copy_block(nocc_*2, nocc_*2+nvirt_, nvirt_*2, nvirtnr_, vfock.slice(nvirt_*2-nvirtnr_, nvirt_*2));
    trans.copy_block(nocc_*2, nocc_*2+nvirt_+nvirtnr_, nvirt_*2, nvirt_-nvirtnr_, vfock.slice(nvirt_, nvirt_*2-nvirtnr_));
  } else {
    if (nclosed_) {
      ZMatrix ofock(ocoeff % fock * ocoeff);
      ofock.diagonalize(eig);
      trans.copy_block(0, 0, nclosed_*2, nclosed_*2, ofock);
    }
    if (nact_ && canonical_) {
      ZMatrix afock = acoeff % fock * acoeff;
      afock.diagonalize(eig);
      trans.copy_block(nclosed_*2, nclosed_*2, nact_*2, nact_*2, afock);
    }
    ZMatrix vfock = vcoeff % fock * vcoeff;
    vfock.diagonalize(eig);
    trans.copy_block(nocc_*2, nocc_*2, nvirt_*2, nvirtnr_*2, vfock.slice((nvirt_-nvirtnr_)*2, nvirt_*2));
    trans.copy_block(nocc_*2, nocc_*2+nvirtnr_*2, nvirt_*2, (nvirt_-nvirtnr_)*2, vfock.slice(0, (nvirt_-nvirtnr_)*2));
  }

  return make_shared<ZCoeff_Block>(*coeff_ * trans, nclosed_, nact_, nvirtnr_, nneg_);
}
