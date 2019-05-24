//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: smith_info.cc
// Copyright (C) 2015 Toru Shiozaki
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
#ifdef COMPILE_SMITH

#include <src/smith/smith_info.h>
#include <src/wfn/zcoeff.h>
#include <src/ci/fci/fci.h>
#include <src/ci/zfci/zharrison.h>

using namespace std;
using namespace bagel;


template<typename DataType>
SMITH_Info<DataType>::SMITH_Info(shared_ptr<const Reference> o, const shared_ptr<const PTree> idata) : ref_(o) {
  stringstream ss;
  method_ = idata->get<string>("method");

  const bool frozen = idata->get<bool>("frozen", true);
  ncore_ = idata->get<int>("ncore", (frozen ? ref_->geom()->num_count_ncore_only()/2 : 0));
  if (ncore_)
    ss << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;
  nfrozenvirt_ = idata->get<int>("nfrozenvirt", 0);
  if (nfrozenvirt_)
    ss << "    * freezing " << nfrozenvirt_ << " orbital" << (nfrozenvirt_^1 ? "s" : "") << " (virtual)" << endl;

  maxiter_ = idata->get<int>("maxiter", 50);
  maxtile_ = idata->get<int>("maxtile", 10);
  cimaxchunk_ = idata->get<int>("cimaxchunk", 317520001);

  do_ms_   = idata->get<bool>("ms",  true);
  do_xms_  = idata->get<bool>("xms", true);
  if (do_xms_ && (method_ == "casa" || method_ == "mrci")) {
    ss << "    * XMS rotation is only appropriate for CASPT2, and will not be used with " << method_ << endl;
    do_xms_ = false;
  }

  sssr_    = idata->get<bool>("sssr", true);
  shift_diag_  = idata->get<bool>("shift_diag", true);
  shift_imag_  = idata->get<bool>("shift_imag", false);
  block_diag_fock_ = idata->get<bool>("block_diag_fock", false);
  orthogonal_basis_ = idata->get<bool>("orthogonal_basis", shift_imag_ ? true : false);

  if (!orthogonal_basis_ && shift_imag_)
    throw runtime_error("Imaginary shift is only applicable when orthogonal basis is used.");

  // check nact() because ciwfn() is nullptr with zero active orbitals
  if (nact() && ciwfn()->nstates() > 1)
    ss << "    * " << (sssr_ ? "SS-SR" : "MS-MR") << " internal contraction is used" << endl;

  // print
  const string sout = ss.str();
  cout << sout << (sout.empty() ? "" : "\n");

  if (ref_->nstate() != 1) {
    if (idata->get<bool>("extract_civectors", false)) {
      vector<int> states = idata->get_vector<int>("extract_state");
      ref_ = extract_ref(states, idata->get<bool>("extract_average_rdms", true));
    } else if (idata->get<bool>("extract_average_rdms", false)) {
      vector<int> states = idata->get_vector<int>("extract_state");
      ref_ = ref_->extract_average_rdm(states);
    }
  }

  // These are not input parameters (set automatically)
  grad_    = idata->get<bool>("_grad", false);

  thresh_ = idata->get<double>("thresh", grad_ ? 1.0e-8 : 1.0e-6);
  shift_  = idata->get<double>("shift", 0.0);
  davidson_subspace_ = idata->get<int>("davidson_subspace", 10);
  thresh_overlap_ = idata->get<double>("thresh_overlap", 1.0e-9);

  // if no convergence is obtained, throw it
  convergence_throw_ = idata->get<bool>("convergence_throw", true);

  // enable restart capability
  restart_ = idata->get<bool>("restart", false);
  restart_each_iter_ = idata->get<bool>("restart_each_iter", restart_);
  state_begin_ = 0;
  restart_iter_ = 0;

  // Restart with MRCI would require us to load amplitudes from previous iterations into DavidsonDiag
  // TODO maybe implement this in the future
  if (restart_ && to_lower(method_) == "mrci")
    throw runtime_error("Restarting is currently only available in SMITH for relativistic perturbation theory methods, not MRCI.");

  // save inputs for pseudospin module
  aniso_data_ = idata->get_child_optional("aniso");
  external_rdm_ = idata->get<string>("external_rdm", "");
  if (nact() && external_rdm_.empty() && !ciwfn()->civectors())
    throw runtime_error("CI vectors are missing. Most likely you ran CASSCF with external RDMs and forgot to specify external_rdm in the smith input block.");
}


template<typename DataType>
SMITH_Info<DataType>::SMITH_Info(shared_ptr<const Reference> o, shared_ptr<const SMITH_Info> info)
  : ref_(o), method_(info->method_), ncore_(info->ncore_), nfrozenvirt_(info->nfrozenvirt_), thresh_(info->thresh_), shift_(info->shift_),
    maxiter_(info->maxiter_), maxtile_(info->maxtile_),
    cimaxchunk_(info->cimaxchunk_), davidson_subspace_(info->davidson_subspace_), grad_(info->grad_),
    do_ms_(info->do_ms_), do_xms_(info->do_xms_), sssr_(info->sssr_),
    shift_diag_(info->shift_diag_), shift_imag_(info->shift_imag_), block_diag_fock_(info->block_diag_fock_), orthogonal_basis_(info->orthogonal_basis_), restart_(info->restart_),
    restart_each_iter_(info->restart_each_iter_), convergence_throw_(info->convergence_throw_), thresh_overlap_(info->thresh_overlap_),
    state_begin_(info->state_begin_), restart_iter_(info->restart_iter_), aniso_data_(info->aniso_data_), external_rdm_(info->external_rdm_) {
}


template<>
tuple<shared_ptr<const RDM<1>>, shared_ptr<const RDM<2>>> SMITH_Info<double>::rdm12(const int ist, const int jst) const {
  FCI_bare fci(ciwfn());
  shared_ptr<const RDM<1>> r1;
  shared_ptr<const RDM<2>> r2;
  if (external_rdm_.empty()) {
    fci.compute_rdm12(ist, jst);
    r1 = fci.rdm1(ist, jst);
    r2 = fci.rdm2(ist, jst);
  } else {
    r1 = fci.read_external_rdm1(ist, jst, external_rdm_);
    r2 = fci.read_external_rdm2(ist, jst, external_rdm_);
  }
  return make_tuple(r1, r2);
}


template<>
tuple<shared_ptr<const RDM<3>>, shared_ptr<const RDM<4>>> SMITH_Info<double>::rdm34(const int ist, const int jst) const {
  FCI_bare fci(ciwfn());
  shared_ptr<const RDM<3>> r3;
  shared_ptr<const RDM<4>> r4;
  if (external_rdm_.empty()) {
    fci.compute_rdm12(ist, jst);
    tie(r3, r4) = fci.rdm34(ist, jst);
  } else {
    r3 = fci.read_external_rdm3(ist, jst, external_rdm_);
    r4 = fci.read_external_rdm4(ist, jst, external_rdm_);
  }
  return make_tuple(r3, r4);
}


template<>
tuple<shared_ptr<const RDM<3>>, shared_ptr<RDM<3>>> SMITH_Info<double>::rdm34f(const int ist, const int jst, shared_ptr<const Matrix> fock) const {
  FCI_bare fci(ciwfn());
  shared_ptr<const RDM<3>> r3;
  shared_ptr<RDM<3>> r4f;
  if (external_rdm_.empty()) {
    fci.compute_rdm12(ist, jst);
    tie(r3, r4f) = fci.rdm34f(ist, jst, fock);
  } else {
    r3 = fci.read_external_rdm3(ist, jst, external_rdm_);
    r4f = fci.read_external_rdm3(ist, jst, external_rdm_, /*fock_contracted=*/true);
  }
  return make_tuple(r3, r4f);
}


template<>
shared_ptr<RDM<3>> SMITH_Info<double>::rdm4f_contract(shared_ptr<const RDM<3>> rdm3, shared_ptr<const RDM<4>> rdm4, shared_ptr<const Matrix> fock) const {
  shared_ptr<RDM<3>> rdm4f = rdm3->clone();

  auto rdm4v = btas::group(group(*rdm4, 6,8), 0,6);
  auto rdm4fv = btas::group(*rdm4f, 0,6);
  contract(1.0, rdm4v, {0,1}, btas::group(*fock,0,2), {1}, 0.0, rdm4fv, {0});

  return rdm4f;
}


template<>
tuple<shared_ptr<const Kramers<2,ZRDM<1>>>, shared_ptr<const Kramers<4,ZRDM<2>>>>
  SMITH_Info<complex<double>>::rdm12(const int ist, const int jst) const {

  ZFCI_bare fci(ciwfn());
  shared_ptr<const Kramers<2,ZRDM<1>>> rdm1;
  shared_ptr<const Kramers<4,ZRDM<2>>> rdm2;
  if (external_rdm_.empty()) {
    rdm1 = fci.rdm1(ist, jst);
    rdm2 = fci.rdm2(ist, jst);
  } else {
    rdm1 = fci.read_external_rdm1(ist, jst, external_rdm_);
    rdm2 = fci.read_external_rdm2(ist, jst, external_rdm_);
  }
  return make_tuple(rdm1, rdm2);
}


template<>
tuple<shared_ptr<const Kramers<6,ZRDM<3>>>, shared_ptr<Kramers<6,ZRDM<3>>>>
  SMITH_Info<complex<double>>::rdm34f(const int ist, const int jst, shared_ptr<const ZMatrix> fock) const {

  ZFCI_bare fci(ciwfn());
  shared_ptr<const Kramers<6,ZRDM<3>>> rdm3;
  shared_ptr<Kramers<6,ZRDM<3>>> rdm4f;

  if (external_rdm_.empty()) {
    tie(rdm3, rdm4f) = fci.rdm34f(ist, jst, fock);
  } else {
    rdm3 = fci.read_external_rdm3(ist, jst, external_rdm_);
    rdm4f = fci.read_external_rdm3(ist, jst, external_rdm_, /*fock_contracted=*/true);
  }

  return make_tuple(rdm3, rdm4f);
}


template<>
shared_ptr<Kramers<6,ZRDM<3>>> SMITH_Info<complex<double>>::rdm4f_contract(shared_ptr<const Kramers<6,ZRDM<3>>> rdm3, shared_ptr<const Kramers<8,ZRDM<4>>> rdm4, shared_ptr<const ZMatrix> fockact) const {
  shared_ptr<Kramers<6,ZRDM<3>>> rdm4f = make_shared<Kramers<6,ZRDM<3>>>();
  const int n = fockact->ndim()/2;

  Kramers<2,ZMatrix> fock;
  fock.emplace(0, fockact->get_submatrix(0, 0, n, n));
  fock.emplace(1, fockact->get_submatrix(n, 0, n, n));
  fock.emplace(2, fockact->get_submatrix(0, n, n, n));
  fock.emplace(3, fockact->get_submatrix(n, n, n, n));
  for (int i = 0; i != 64; ++i) {
    for (int j = 0; j != 4; ++j) {
      auto work = make_shared<ZRDM<3>>(n);
      shared_ptr<const ZMatrix> cfock = fock.at(j);
      shared_ptr<const ZRDM<4>> crdm = rdm4->get_data(i * 4 + j);
      if (!crdm) continue;

      auto wgr = btas::group(*work, 0,6);
      auto crdmgr = btas::group(btas::group(*crdm, 6,8),0,6);
      auto fgr = btas::group(*cfock, 0,2);
      btas::contract(1.0, crdmgr, {0,1}, fgr, {1}, 0.0, wgr, {0});
      rdm4f->add(i, work);
    }
  }

  return rdm4f;
}


template<>
tuple<shared_ptr<const Kramers<6,ZRDM<3>>>, shared_ptr<const Kramers<8,ZRDM<4>>>>
  SMITH_Info<complex<double>>::rdm34(const int ist, const int jst) const {

  ZFCI_bare fci(ciwfn());
  shared_ptr<const Kramers<6,ZRDM<3>>> rdm3;
  shared_ptr<const Kramers<8,ZRDM<4>>> rdm4;
  if (external_rdm_.empty()) {
    rdm3 = fci.rdm3(ist, jst);
    rdm4 = fci.rdm4(ist, jst);
  } else {
    rdm3 = fci.read_external_rdm3(ist, jst, external_rdm_);
    rdm4 = fci.read_external_rdm4(ist, jst, external_rdm_);
  }
  return make_tuple(rdm3, rdm4);
}


template<>
shared_ptr<const RDM<1>> SMITH_Info<double>::rdm1_av() const {
  return ref_->rdm1_av();
}


template<>
shared_ptr<const ZRDM<1>> SMITH_Info<complex<double>>::rdm1_av() const {
  return nullptr;
}


template<>
shared_ptr<const CIWfn> SMITH_Info<double>::ciwfn() const {
  return ref_->ciwfn();
}


template<>
shared_ptr<const RelCIWfn> SMITH_Info<complex<double>>::ciwfn() const {
  return dynamic_pointer_cast<const RelReference>(ref_)->ciwfn();
}


template<>
shared_ptr<const Matrix> SMITH_Info<double>::coeff() const {
  return ref_->coeff();
}


template<>
shared_ptr<const ZMatrix> SMITH_Info<complex<double>>::coeff() const {
  shared_ptr<const ZCoeff_Striped> c = dynamic_pointer_cast<const RelReference>(ref_)->relcoeff();
  return c->block_format(nclosed(), nact(), nvirt()+nfrozenvirt(), 0);
}


template<>
shared_ptr<const Matrix> SMITH_Info<double>::hcore() const {
  return ref_->hcore();
}


template<>
shared_ptr<const ZMatrix> SMITH_Info<complex<double>>::hcore() const {
  // TODO implement
  assert(false);
  return nullptr;
}


template<typename DataType>
shared_ptr<const Reference>  SMITH_Info<DataType>::extract_ref(const vector<int> states, const bool extract_rdm) const {
  shared_ptr<const Reference> out = ref_;
  cout << "    * Running " << (do_xms_ ? "X" : "" ) << (do_ms_ ? "MS " : "") << method_ << " for all retained states from a multi-state reference." << endl;
  out = ref_->extract_state(states, extract_rdm);
  return out;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class bagel::SMITH_Info<double>;
template class bagel::SMITH_Info<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_CLASS_EXPORT_IMPLEMENT(SMITH_Info<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(SMITH_Info<complex<double>>)

#endif
