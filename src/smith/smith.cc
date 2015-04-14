//
// BAGEL - Parallel electron correlation program.
// Filename: smith.cc
// Copyright (C) 2013 Matthew MacLeod
//
// Author: Matthew K. MacLeod <matthew.macleod@northwestern.edu>
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


#include <bagel_config.h>
#include <src/smith/smith.h>
#include <src/smith/MRCI.h>
#include <src/smith/CASPT2.h>


using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

Smith::Smith(const shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : Method(idata, g, r) {
  const string method = to_lower(idata_->get<string>("method", "caspt2"));

  // make a smith_info class
  auto info = make_shared<SMITH_Info<double>>(r, idata);

#ifdef COMPILE_SMITH
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
  algo_->solve();

#ifdef COMPILE_SMITH
  if (algo_->info()->grad()) {
    auto algop = dynamic_pointer_cast<CASPT2::CASPT2>(algo_);
    assert(algop);

    algop->solve_deriv();
    dm1_ = algop->rdm12();
    dm11_ = algop->rdm11();
    dm2_ = algop->rdm21();

    // compute <1|1>
    wf1norm_ = algop->correlated_norm();

    // convert ci derivative tensor to civec
    cider_ = algop->ci_deriv(ref_->ciwfn()->det());

    // todo check
    coeff_ = make_shared<Coeff>(*algop->coeff());
  }
#endif
}
