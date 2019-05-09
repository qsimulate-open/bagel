//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: smith.cc
// Copyright (C) 2013 Toru Shiozaki
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
  // print a header
  if (!idata->get<bool>("_grad", false))
    cout << "  === SMITH program ===" << endl << endl;

  // make a smith_info class
  auto info = make_shared<const SMITH_Info<double>>(r, idata);

  if (idata->get<bool>("_grad", false) && !info->do_ms() && r->nstate() > 1)
    throw runtime_error("CASPT2 nuclear gradients are not implemented for SS-CASPT2 with multiple reference states");

  {
    // print memory requirements
    const size_t nact = r->nact();
    const size_t nstate = r->nstate();
    const size_t nclosed = r->nclosed() - info->ncore();
    const size_t nvirtual = r->nvirt() - info->nfrozenvirt();
    const size_t nocc = nclosed + nact;
    const size_t norb = nocc + nvirtual;
    const size_t rsize = nclosed*nocc*nvirtual*nvirtual + nclosed*nclosed*nact*(nvirtual+nact) + pow(nact,2)*(norb*nvirtual+nact*nclosed);
    cout << "    * Approximate memory requirement for SMITH calculation per MPI process:" << endl;
    cout << "      o Storage requirement for T-amplitude, lambda, and residual is ";
    if (info->orthogonal_basis()) {
      // During iteration      : davidson_subspace * (Lambda or T (ortho) + Residual (ortho)) * nstate + (T-amp(orthogonal, redundant) + Res + Lambda) * nstate * nstate
      cout << setprecision(2) << (info->davidson_subspace() * 2 + nstate * (4 + info->grad() ? 2 : 0)) * rsize * (info->sssr() ? 1 : nstate) * 8.e-9 / mpi__->size() << " GB" << endl;
    } else {
      // During iteration      : davidson_subspace * (Lambda or T (redun) + Residual (r)) * nstate + (T-amp(r) + Res + Lambda) * nstate * nstate
      cout << setprecision(2) << (info->davidson_subspace() * 2 + nstate * (2 + info->grad() ? 1 : 0)) * rsize * (info->sssr() ? 1 : nstate) * 8.e-9 / mpi__->size() << " GB" << endl;
    }
    cout << "      o Storage requirement for MO integrals is ";
    cout << setprecision(2) << (norb*norb*2 + nocc*nocc*(nact+nvirtual)*(nact+nvirtual)) * 8.e-9 / mpi__->size() << " GB" << endl;
    if (info->grad()) {
      cout << "      o Storage requirement for SMITH-computed gradient tensors is ";
      cout << setprecision(2) << nstate*nstate*(pow(nact,6)*2 + pow(nact,4) + pow(nact,2) + 1) * 8.e-9 << " GB" << endl;
    }
  }

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
#else
  throw logic_error("You must enable SMITH during compilation for this method to be available.");
#endif
}


void Smith::compute_gradient(const int istate, const int jstate, shared_ptr<const NacmType> nacmtype, const bool nocider) {
#ifdef COMPILE_SMITH
  if (!(algo_->info()->grad())) {
    auto algop = dynamic_pointer_cast<CASPT2::CASPT2>(algo_);
    algop->solve_dm(istate, jstate);
    msrot_ = algop->msrot();
    coeff_ = algop->coeff();
    vd1_ = algop->vden1();
  } else {
    auto algop = make_shared<CASPT2::CASPT2>(*(dynamic_pointer_cast<CASPT2::CASPT2>(algo_)));
    assert(algop);

    algop->solve_gradient(istate, jstate, nacmtype, nocider);
    dm1_ = algop->rdm12();
    dm11_ = algop->rdm11();
    dm2_ = algop->rdm21();
    dcheck_ = algop->dcheck();
    if (istate != jstate)
      vd1_ = algop->vden1();

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
  }
#else
  throw logic_error("You must enable SMITH during compilation for this method to be available.");
#endif
}


RelSmith::RelSmith(const shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : Method(idata, g, r) {
#ifdef COMPILE_SMITH
  // print a header
  if (!idata->get<bool>("_grad", false))
    cout << "  === SMITH program ===" << endl << endl;

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
#ifndef DISABLE_SERIALIZATION
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
#else
  throw logic_error("You must enable serialization to make \"continue\" to be available.");
#endif
  }
#else
  throw logic_error("You must enable SMITH during compilation for this method to be available.");
#endif
}
