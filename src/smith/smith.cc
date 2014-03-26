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


#include <src/smith/smith.h>
#include <src/smith/MP2.h>
#include <src/smith/CAS_all_active.h>
#include <src/smith/CAS_test.h>


using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

Smith::Smith(const shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : Method(idata, g, r) {
  string method = idata_->get<string>("method", "mp2");

  // make a smith_info class
  auto info = make_shared<SMITH_Info>(r, idata);

  if (method == "mp2") {
    algo_ = make_shared<MP2::MP2<Storage_Incore>>(info);
  } else if (method == "caspt2") {
    algo_ = make_shared<CAS_all_active::CAS_all_active<Storage_Incore>>(info);
  } else if (method == "caspt2-test") {
    algo_ = make_shared<CAS_test::CAS_test<Storage_Incore>>(info);
  } else {
    stringstream ss; ss << method << " method is not implemented in SMITH";
    throw logic_error(ss.str());
  }
}

void Smith::compute() {
  algo_->solve();

  dm1_ = dynamic_pointer_cast<CAS_test::CAS_test<Storage_Incore>>(algo_)->rdm1();
  dm2_ = dynamic_pointer_cast<CAS_test::CAS_test<Storage_Incore>>(algo_)->rdm2();

  // calculate unrelaxed dipole moment from correlated dm
  correction_ = dynamic_pointer_cast<CAS_test::CAS_test<Storage_Incore>>(algo_)->rdm1_correction();
  algo_->dipole(dm1_,correction_,"CASPT2 Unrelaxed").compute();

  // convert ci derivative tensor to civec
  cider_ = dynamic_pointer_cast<CAS_test::CAS_test<Storage_Incore>>(algo_)->ci_deriv();

  // todo check
  coeff_ = dynamic_pointer_cast<CAS_test::CAS_test<Storage_Incore>>(algo_)->coeff();

  cout << "  * Printing ci derivative civec:" << endl;
  cider_->print(0.1e-15);
  cout << "  * Printing civec ci derivative * cI =     " <<  setprecision(10) << cider_->dot_product(*(algo_->rdm0deriv())) << endl;


}
