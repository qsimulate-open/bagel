//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbital/asd_orbopt.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
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


#include <fstream>
#include <src/scf/hf/fock.h>
#include <src/asd/orbital/asd_orbopt.h>

#include <src/wfn/construct_method.h>

using namespace std;
using namespace bagel;

ASD_OrbOpt::ASD_OrbOpt(shared_ptr<const PTree> idat, shared_ptr<Dimer> dimer)
  : Method(idat, dimer->sgeom(), dimer->sref()), dimer_(dimer), hcore_(make_shared<Hcore>(dimer->sgeom())) {
  common_init();
}


void ASD_OrbOpt::common_init() {
  print_header();

  // first set coefficient
  coeff_ = ref_->coeff();

  // get maxiter from the input
  max_iter_ = idata_->get<int>("maxiter", 50);
  // get nstate from the ASD input
  nstate_ = idata_->get_child_optional("asd")->get<int>("nstates", 1);
  // get thresholds
  gradient_thresh_ = idata_->get<double>("gradient_thresh", 1.0e-4);
  rotation_thresh_ = idata_->get<double>("rotation_thresh", 1.0e-4);
  energy_thresh_ = idata_->get<double>("energy_thresh", 1.0e-6);

  //semi_canonical orbital at convergence
  semi_canonical_ = idata_->get<bool>("semi_canonical", false);

  // active space
  nact_ = ref_->nact();
  nactA_ = dimer_->active_refs().first->nact();
  nactB_ = dimer_->active_refs().second->nact();
  nactcloA_ = dimer_->isolated_refs().first->nclosed() - dimer_->active_refs().first->nclosed();
  nactcloB_ = dimer_->isolated_refs().second->nclosed() - dimer_->active_refs().second->nclosed();
  rasA_ = {0, nactA_, 0};
  rasB_ = {0, nactB_, 0};
  if (idata_->get_child("asd")->get<string>("method") == "ras") {
    auto restrictions = idata_->get_child("asd")->get_child("restricted");
    auto get_restricted_data = [] (shared_ptr<const PTree> i) { return i->get_array<int, 3>("orbitals"); };

    if (restrictions->size() == 1) {
      rasA_ = get_restricted_data(*restrictions->begin());
      rasB_ = rasA_;
    } else if (restrictions->size() == 2) {
      auto iter = restrictions->begin();
      rasA_ = get_restricted_data(*iter++);
      rasB_ = get_restricted_data(*iter);
    } else {
      throw logic_error("One or two sets of restrictions must be provided.");
    }
  }

  assert(nactA_ + nactB_ == nact_);

  nclosed_ = ref_->nclosed();

  nocc_ = nclosed_ + nact_;

  nbasis_ = coeff_->mdim();
  nvirt_ = nbasis_ - nocc_;
  if (nvirt_ < 0) throw runtime_error("It appears that nvirt < 0. Check the nocc value");

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    *  unit A  : " << setw(6) << nactA_ << " (" << rasA_[0] << "," << rasA_[1] << "," << rasA_[2] << ")" << endl;
  cout << "    *  unit B  : " << setw(6) << nactB_ << " (" << rasB_[0] << "," << rasB_[1] << "," << rasB_[2] << ")" << endl;
  cout << "    * nocc     : " << setw(6) << nocc_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl;

  const int idel = geom_->nbasis() - nbasis_;
  if (idel)
    cout << "      Due to linear dependency, " << idel << (idel==1 ? " function is" : " functions are") << " omitted" << endl;

  cout <<  "  === ASD Orbital Optimization iteration (" + geom_->basisfile() + ") ===" << endl << endl;

}


ASD_OrbOpt::~ASD_OrbOpt() {

}


void ASD_OrbOpt::print_header() const {
  cout << "  --------------------------------------------" << endl;
  cout << "      ASD Orbital Optimization calculation    " << endl;
  cout << "  --------------------------------------------" << endl << endl;
}

void ASD_OrbOpt::print_iteration(int iter, int miter, int tcount, const vector<double> energy, const double error, double max_r, double delta_e, const double time, const int itype, bool fix_ci) const {
  auto print_iteration_type = [] (const int i) {
    string out;
    if (i == 0) out = "intra";
    else if (i == 1) out = "inter";
    else out = "full";
    return out;
  };

  if (energy.size() != 1 && iter) cout << endl;

  if (iter == 0) cout << setw(7) << "iter" << setw(31) << "energy" << setw(13) << "gradient" << setw(10) << "rotation" << setw(10) << "delta_e" << endl;

  int i = 0;
  for (auto& e : energy) {
    cout << "  " << setw(5) << iter << setw(3) << i << setw(4) << miter << setw(4) << tcount
                 << setw(20) << fixed << setprecision(12) << e << "   "
                 << setw(10) << scientific << setprecision(2) << (i==0 ? error : 0.0)
                 << setw(10) << scientific << setprecision(2) << (i==0 ? max_r : 0.0)
                 << setw(10) << scientific << setprecision(2) << (i==0 ? delta_e : 0.0)
                 << fixed << setw(10) << setprecision(2) << time
                 << setw(10) << print_iteration_type(itype)
                 << setw(10) << (fix_ci ? "fix" : " ") << endl;
    ++i;
  }
}


static streambuf* backup_stream_;
static ofstream* ofs_;


void ASD_OrbOpt::mute_stdcout() {
  ofstream* ofs(new ofstream("asd_orbopt.log",(backup_stream_ ? ios::app : ios::trunc)));
  ofs_ = ofs;
  backup_stream_ = cout.rdbuf(ofs->rdbuf());
}


void ASD_OrbOpt::resume_stdcout() {
  cout.rdbuf(backup_stream_);
  delete ofs_;
}


shared_ptr<Matrix> ASD_OrbOpt::Qvec(const int n, const int m, shared_ptr<const Matrix> coeff, const size_t nclosed) const {
  assert(n == coeff->mdim());

  shared_ptr<DFHalfDist> half;
  const MatView cdata = coeff_->slice(nclosed_, nclosed_+nact_);
  half = geom_->df()->compute_half_transform(cdata);

  // J^{-1}(D|xy)
  // TODO : DFDistT needs to be modified to handle cases where number of nodes is larger than half->nocc() * cdata.mdim()
  shared_ptr<const DFFullDist> full;
  if (half->nocc() * coeff->mdim() > mpi__->size()) {
    full = half->apply_JJ()->compute_second_transform(coeff->slice(nclosed, nclosed+m));
  } else {
    full = half->compute_second_transform(coeff->slice(nclosed, nclosed+m))->apply_JJ();
  }

  // [D|tu] = (D|xy)Gamma_xy,tu
  shared_ptr<const DFFullDist> prdm = full->apply_2rdm(*rdm2_);

  // (r,u) = (rt|D)[D|tu]
  shared_ptr<const Matrix> tmp = half->form_2index(prdm, 1.0);

  // MO transformation of the first index
  auto out = make_shared<Matrix>(*coeff % *tmp);
  assert(n == out->ndim() && m == out->mdim());
  return out;

}


double ASD_OrbOpt::check_symmetric(shared_ptr<Matrix>& mat) const {
  int n = mat->ndim();
  assert(n == mat->mdim());
  auto tran = make_shared<Matrix>(n,n);
  tran = mat->transpose();
  auto subt = make_shared<Matrix>(n,n);
  *subt = *mat - *tran;
  return subt->rms();
}


shared_ptr<const Coeff> ASD_OrbOpt::update_coeff(const shared_ptr<const Matrix> cold, shared_ptr<const Matrix> mat) const {
  auto cnew = make_shared<Coeff>(*cold);
  int nbas = cold->ndim();
  assert(nbas == geom_->nbasis());
  dgemm_("N", "N", nbas, nact_, nact_, 1.0, cold->data()+nbas*nclosed_, nbas, mat->data(), nact_,
                   0.0, cnew->data()+nbas*nclosed_, nbas);
  return cnew;
}


shared_ptr<const Reference> ASD_OrbOpt::conv_to_ref() const {
  cout << "ASD_OrbOpt: conv_to_ref: not yet implemented" << endl;
  //TODO
  assert(false);
  return nullptr;
}


shared_ptr<const Coeff> ASD_OrbOpt::semi_canonical_orb() const {
  auto rdm1mat = make_shared<Matrix>(*rdm1_->rdm1_mat(/*nclosed*/0));

  auto ocoeff = coeff_->slice_copy(0, nclosed_);
  auto acoeff = coeff_->slice_copy(nclosed_, nocc_);

  shared_ptr<Matrix> semi_coeff = coeff_->copy();
  //ncloesd_ nact_ nvirt_ nbasis_, nactA_, nactB_, nactcloA_, nactcloB_
  cout << "Semi-canonical transformation" << endl;
  cout << "nclosed  : " << nclosed_ << endl;
  cout << "nact     : " << nact_    << endl;
  cout << "nvirt    : " << nvirt_   << endl;
  cout << "nbasis   : " << nbasis_  << endl;
  cout << "nactA    : " << nactA_   << endl;
  cout << "nactB    : " << nactB_   << endl;
  cout << "nactcloA    : " << nactcloA_   << endl;
  cout << "nactcloB    : " << nactcloB_   << endl;
  {//Monomer A
    //core
    auto acoeff_A = acoeff->slice_copy(0,nactA_);
    auto ccoeff_A = make_shared<Matrix>(nbasis_, nclosed_+nactcloB_); //nclosed : shared closed including closed activeB
    ccoeff_A->copy_block(0,0, nbasis_,nclosed_, ocoeff);
    ccoeff_A->copy_block(0,nclosed_, nbasis_,nactcloB_, acoeff->get_submatrix(0,nactA_, nbasis_,nactcloB_)); //embed activeB
    shared_ptr<const Matrix> ofockao_A = make_shared<Fock<1>>(geom_, ref_->hcore(), nullptr, ccoeff_A, /*store*/false, /*rhf*/true);
    //active
    auto rdm1_A = make_shared<Matrix>(*rdm1mat->get_submatrix(0,0, nactA_,nactA_));
    rdm1_A->sqrt();
    auto acoeffw_A = make_shared<Matrix>(*acoeff_A * (1.0/sqrt(2.0)) * *rdm1_A);
    auto fockao_A = make_shared<Fock<1>>(geom_, ofockao_A, nullptr, acoeffw_A, /*store*/false, /*rhf*/true);
    // MO Fock
    VectorB eigs(nactA_);
    auto fockact = make_shared<Matrix>(*acoeff_A % *fockao_A  * *acoeff_A);
    fockact->diagonalize(eigs);
    for (int i = 0; i < nactA_; ++i)
      cout << i << "(A) = " << eigs[i] << endl;
    *acoeff_A *= *fockact;

    size_t A_position = nclosed_; //for A
    for (int i = 0; i < nactA_; ++i)
      copy_n(acoeff_A->element_ptr(0, i), nbasis_, semi_coeff->element_ptr(0,A_position++));
  }
  {//Monomer B
    //core
    auto acoeff_B = acoeff->slice_copy(nactA_,nact_);
    auto ccoeff_B = make_shared<Matrix>(nbasis_, nclosed_+nactcloA_);
    ccoeff_B->copy_block(0,0, nbasis_,nclosed_, ocoeff);
    ccoeff_B->copy_block(0,nclosed_, nbasis_,nactcloA_, acoeff->get_submatrix(0,0, nbasis_,nactcloA_));
    shared_ptr<const Matrix> ofockao_B = make_shared<Fock<1>>(geom_, ref_->hcore(), nullptr, ccoeff_B, /*store*/false, /*rhf*/true);
    //active
    auto rdm1_B = make_shared<Matrix>(*rdm1mat->get_submatrix(nactA_,nactA_, nactB_,nactB_));
    rdm1_B->sqrt();
    auto acoeffw_B = make_shared<Matrix>(*acoeff_B * (1.0/sqrt(2.0)) * *rdm1_B);
    auto fockao_B = make_shared<Fock<1>>(geom_, ofockao_B, nullptr, acoeffw_B, /*store*/false, /*rhf*/true);
    // MO Fock
    VectorB eigs(nactB_);
    auto fockact = make_shared<Matrix>(*acoeff_B % *fockao_B  * *acoeff_B);
    fockact->diagonalize(eigs);
    for (int i = 0; i < nactB_; ++i)
      cout << i << "(B) = " << eigs[i] << endl;
    *acoeff_B *= *fockact;

    size_t B_position = nclosed_+nactA_;
    for (int i = 0; i < nactB_; ++i)
      copy_n(acoeff_B->element_ptr(0, i), nbasis_, semi_coeff->element_ptr(0,B_position++));
  }

  return make_shared<Coeff>(*semi_coeff);
}
