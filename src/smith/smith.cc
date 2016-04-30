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
using namespace bagel::SMITH;
#endif
using namespace std;
using namespace bagel;

Smith::Smith(const shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : Method(idata, g, r) {
  const string method = to_lower(idata_->get<string>("method", "caspt2"));

#ifdef COMPILE_SMITH
  // make a smith_info class
  auto info = make_shared<SMITH_Info<double>>(r, idata);

  if (method == "caspt2") {
    algo_ = make_shared<CASPT2::CASPT2>(info);
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
    uwumat_ = algop->uwumat();

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
      sdm1_ = make_shared<Matrix>(*sp->rdm12() * 2.0 - *dm1_);
      sdm11_ = make_shared<Matrix>(*sp->rdm11() * 2.0 - *dm11_);
    }
  }
#endif
}

RelSmith::RelSmith(const shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : Method(idata, g, r) {
  const string method = to_lower(idata_->get<string>("method", "caspt2"));

#ifdef COMPILE_SMITH
  // make a smith_info class
  auto info = make_shared<SMITH_Info<complex<double>>>(r, idata);

  if (method == "caspt2") {
    algo_ = make_shared<RelCASPT2::RelCASPT2>(info);
  } else if (method == "mrci") {
    algo_ = make_shared<RelMRCI::RelMRCI>(info);
  } else {
    stringstream ss; ss << method << " method is not implemented in SMITH";
    throw logic_error(ss.str());
  }
#endif
}
