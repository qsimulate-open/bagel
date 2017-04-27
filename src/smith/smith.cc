//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: smith.cc
// Copyright (C) 2013 Matthew MacLeod
//
// Author: Matthew K. MacLeod <matthew.macleod@northwestern.edu>
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

#include <src/smith/smith.h>

#ifdef COMPILE_SMITH
#include <src/smith/mrci/MRCI.h>
#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/SPCASPT2.h>
#include <src/smith/relcaspt2/RelCASPT2.h>
#include <src/smith/casa/CASA.h>
#include <src/smith/relcasa/RelCASA.h>
using namespace bagel::SMITH;
#endif
using namespace std;
using namespace bagel;

Smith::Smith(const shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : Method(idata, g, r) {
  const string method = to_lower(idata_->get<string>("method", "caspt2"));

#ifdef COMPILE_SMITH
  // make a smith_info class
  auto info = make_shared<const SMITH_Info<double>>(r, idata);
  if (info->restart())
    cout << " ** Restarting calculations is currently unavailable for non-relativistic SMITH methods.  Serialized archives will not be generated." << endl << endl;

  if (method == "caspt2") {
    algo_ = make_shared<CASPT2::CASPT2>(info);
  } else if (method == "casa") {
    algo_ = make_shared<CASA::CASA>(info);
  } else if (method == "mrci") {
    algo_ = make_shared<MRCI::MRCI>(info);
  } else {
#else
  {
#endif
    stringstream ss; ss << method << " method is not implemented in SMITH";
    throw logic_error(ss.str());
  }
}

void Smith::compute() {
#ifdef COMPILE_SMITH
  algo_->solve();

  if (algo_->info()->grad()) {
    auto algop = dynamic_pointer_cast<CASPT2::CASPT2>(algo_);
    assert(algop);

    algop->solve_deriv();
    dm1_ = algop->rdm12();
    dm11_ = algop->rdm11();
    dm2_ = algop->rdm21();
    dcheck_ = algop->dcheck();

    // compute <1|1>
    wf1norm_ = algop->correlated_norm();
    // convert ci derivative tensor to civec
    cider_ = algop->ci_deriv();
    msrot_ = algop->msrot();
    coeff_ = algop->coeff();

    // if spin-density is requested...
    if (idata_->get<bool>("_hyperfine")) {
      auto sp = make_shared<SPCASPT2::SPCASPT2>(*algop);
      sp->solve();
      sdm1_ = make_shared<Matrix>(*sp->rdm12() * 2.0 - *dm1_); // CAUTION! dm1 includes <1|1>D0 where as sp->rdm12() does not
      sdm11_ = make_shared<Matrix>(*sp->rdm11() * 2.0 - *dm11_);
    }
  } else if (algo_->info()->nacm()) {
    auto algop = dynamic_pointer_cast<CASPT2::CASPT2>(algo_);
    assert(algop);

    algop->solve_nacme();
    dm1_ = algop->rdm12();
    dm11_ = algop->rdm11();
    vd1_ = algop->vden1();
    dm2_ = algop->rdm21();
    dcheck_ = algop->dcheck();

    // compute <1|1>
    wf1norm_ = algop->correlated_norm();
    // convert ci derivative tensor to civec
    cider_ = algop->ci_deriv();
    foeig_ = algop->e0all();
    xmsrot_ = algop->xmsrot();
    heffrot_ = algop->heffrot();
    msrot_ = algop->msrot();
    coeff_ = algop->coeff();
  } else if (algo_->info()->method() == "caspt2") {
    auto algop = dynamic_pointer_cast<CASPT2::CASPT2>(algo_);
    if (algo_->info()->target2() != -1) {
      algop->solve_dm();
      msrot_ = algop->msrot();
      xmsrot_ = algop->xmsrot();
      heffrot_ = algop->heffrot();
      coeff_ = algop->coeff();
      vd1_ = algop->vden1();
    }
  }
#else
  throw logic_error("You must enable SMITH during compilation for this method to be available.");
#endif
}

RelSmith::RelSmith(const shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : Method(idata, g, r) {
#ifdef COMPILE_SMITH
  const string method = to_lower(idata_->get<string>("method", "caspt2"));
  if (!dynamic_pointer_cast<const RelReference>(r) && method != "continue")
    throw runtime_error("Relativistic correlation methods require a fully relativistic reference wavefunction.");

  if (method != "continue") {
    // make a smith_info class
    auto info = make_shared<const SMITH_Info<complex<double>>>(r, idata);
    if (method == "caspt2") {
      algo_ = make_shared<RelCASPT2::RelCASPT2>(info);
    } else if (method == "casa") {
      algo_ = make_shared<RelCASA::RelCASA>(info);
    } else if (method == "mrci") {
      algo_ = make_shared<RelMRCI::RelMRCI>(info);
    } else {
      stringstream ss; ss << method << " method is not implemented in RelSMITH";
      throw logic_error(ss.str());
    }

  } else {
    // method == "continue" - so load SMITH_Info and T2 amplitudes from Archives
    Timer mtimer;

    shared_ptr<const SMITH_Info<complex<double>>> info;
    string arch_prefix = idata_->get<string>("archive_location", "");
    if (arch_prefix.size() > 2 && arch_prefix.back() != '/')
      arch_prefix += "/";
    {
      IArchive archive(arch_prefix + "RelSMITH_info");
      shared_ptr<SMITH_Info<complex<double>>> ptr;
      archive >> ptr;
      const int state_begin = idata_->get<int>("state_begin", 0);
      const int restart_iter = idata_->get<int>("restart_iter", 0);
      ptr->set_restart_params(state_begin, restart_iter);
      info = shared_ptr<SMITH_Info<complex<double>>>(ptr);
    }
    mtimer.tick_print("Load RelSMITH info Archive");

    ref_ = info->ref();
    geom_ = ref_->geom();
    const string method = to_lower(info->method());
    assert(method == "caspt2" || method == "casa");
    if (method == "caspt2") {
      arch_prefix += "RelCASPT2";
    } else if (method == "casa") {
      arch_prefix += "RelCASA";
    }

    if (method == "caspt2")
      algo_ = make_shared<RelCASPT2::RelCASPT2>(info);
    else if (method == "casa")
      algo_ = make_shared<RelCASA::RelCASA>(info);
    mtimer.tick_print("Construct " + method + " architecture");

    for (int ist = 0; ist <= info->state_begin(); ++ist) {
      if (ist < info->state_begin() || info->restart_iter() > 0) {
        string arch = arch_prefix + "_t2_" + to_string(ist);
        if (ist < info->state_begin())
           arch += "_converged";
        else
           arch += "_iter_" + to_string(info->restart_iter());
        IArchive archive(arch);
        shared_ptr<MultiTensor_<complex<double>>> t2in;
        archive >> t2in;
        if (method == "caspt2")
          (dynamic_pointer_cast<RelCASPT2::RelCASPT2>(algo_))->load_t2all(t2in, ist);
        else
          (dynamic_pointer_cast<RelCASA::RelCASA>(algo_))->load_t2all(t2in, ist);
      }
    }
    mtimer.tick_print("Load T-amplitude Archive (RelSMITH)");
  }
#else
  throw logic_error("You must enable SMITH during compilation for this method to be available.");
#endif
}
