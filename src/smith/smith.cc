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
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

Smith::Smith(const boost::property_tree::ptree& idata, shared_ptr<const Reference> ref) {
  string method = idata.get<string>("method", "mp2");
  if (method == "mp2") {
    algo_ = make_shared<MP2::MP2<Storage_Incore>>(ref);
  } else if (method == "caspt2") {
    algo_ = make_shared<CAS_all_active::CAS_all_active<Storage_Incore>>(ref);
  } else if (method == "caspt2-test") {
    algo_ = make_shared<CAS_test::CAS_test<Storage_Incore>>(ref);
  } else {
    stringstream ss; ss << method << " method is not implemented in SMITH";
    throw logic_error(ss.str());
  }
}

void Smith::compute() {
  algo_->solve();

  shared_ptr<const Matrix> dm1 = algo_->rdm1();
  dm1->print("dm1", 20);
  cout << dm1->ndim() << " " << dm1->mdim() << endl;

  // calculate unrelaxed dipole moment from dm
  algo_->dipole().compute();

}
