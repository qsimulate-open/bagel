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
#include <src/wfn/relcoeff.h>
#include <src/ci/fci/fci.h>
#include <src/ci/zfci/zharrison.h>

using namespace std;
using namespace bagel;


template<typename DataType>
SMITH_Info<DataType>::SMITH_Info(shared_ptr<const Reference> o, const shared_ptr<const PTree> idata) : ref_(o) {
  method_ = idata->get<string>("method");

  const bool frozen = idata->get<bool>("frozen", true);
  ncore_ = idata->get<int>("ncore", (frozen ? ref_->geom()->num_count_ncore_only()/2 : 0));
  if (ncore_)
    cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;
  nfrozenvirt_ = idata->get<int>("nfrozenvirt", 0);
  if (nfrozenvirt_)
    cout << "    * freezing " << nfrozenvirt_ << " orbital" << (nfrozenvirt_^1 ? "s" : "") << " (virtual)" << endl;

  maxiter_ = idata->get<int>("maxiter", 50);
  maxtile_ = idata->get<int>("maxtile", 10);
  cimaxtile_ = idata->get<int>("cimaxtile", 10);

  do_ms_   = idata->get<bool>("ms",  true);
  do_xms_  = idata->get<bool>("xms", false);
  if (do_xms_ && (method_ == "casa" || method_ == "mrci")) {
    cout << "    * XMS rotation is only appropriate for CASPT2, and will not be used with " << method_ << endl;
    do_xms_ = false;
  }

  sssr_    = idata->get<bool>("sssr", false);
  shift_diag_  = idata->get<bool>("shift_diag", true);
  block_diag_fock_ = idata->get<bool>("block_diag_fock", false);

  // check nact() because ciwfn() is nullptr with zero active orbitals
  if (nact() && ciwfn()->nstates() > 1)
    cout << "    * " << (sssr_ ? "SS-SR" : "MS-MR") << " internal contraction is used" << endl;

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
  target_  = idata->get<int>("_target", -1);
  grad_    = idata->get<bool>("_grad", false);
  target2_ = idata->get<int>("_target2", -1);
  nacm_    = idata->get<bool>("_nacm", false);
  nacmtype_= idata->get<int>("_nacmtype", 0);

  thresh_ = idata->get<double>("thresh", (grad_||nacm_) ? 1.0e-8 : 1.0e-6);
  shift_  = idata->get<double>("shift", 0.0);
  davidson_subspace_ = idata->get<int>("davidson_subspace", 10);
  thresh_overlap_ = idata->get<double>("thresh_overlap", 1.0e-9);

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
  if (external_rdm_.empty() && !ciwfn()->civectors())
    throw runtime_error("CI vectors are missing. Most likely you ran CASSCF with external RDMs and forgot to specify external_rdm in the smith input block.");  
  if (!external_rdm_.empty() && is_same<DataType,double>::value) 
    throw logic_error("so far the external RDMs are only interfaced to relativistic theories. TODO");

  assert(!(grad_ && target_ < 0));
  assert(!(nacm_ && target2_ < 0));
}


template<typename DataType>
SMITH_Info<DataType>::SMITH_Info(shared_ptr<const Reference> o, shared_ptr<const SMITH_Info> info)
  : ref_(o), method_(info->method_), ncore_(info->ncore_), nfrozenvirt_(info->nfrozenvirt_), thresh_(info->thresh_), shift_(info->shift_),
    maxiter_(info->maxiter_), target_(info->target_), target2_(info->target2_), nacmtype_(info->nacmtype_),
    maxtile_(info->maxtile_), cimaxtile_(info->cimaxtile_), davidson_subspace_(info->davidson_subspace_), grad_(info->grad_), nacm_(info->nacm_),
    do_ms_(info->do_ms_), do_xms_(info->do_xms_), sssr_(info->sssr_),
    shift_diag_(info->shift_diag_), thresh_overlap_(info->thresh_overlap_), aniso_data_(info->aniso_data_) {
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
  shared_ptr<const RelCoeff_Striped> c = dynamic_pointer_cast<const RelReference>(ref_)->relcoeff();
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
template class SMITH_Info<double>;
template class SMITH_Info<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_CLASS_EXPORT_IMPLEMENT(SMITH_Info<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(SMITH_Info<complex<double>>)

#endif
