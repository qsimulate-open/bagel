//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hess.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Bess Vlaisavljevich <bess.vlaisavljevich@northwestern.edu>
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

#include <string>
#include <src/grad/hess.h>
#include <src/wfn/construct_method.h>
#include <src/grad/force.h>
#include <src/util/math/xyzfile.h>

using namespace std;
using namespace bagel;

Hess::Hess(shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : idata_(idata), geom_(g), ref_(r) {

}

void Hess::compute() {
  auto input = idata_->get_child("method");
  const string jobtitle = to_lower(idata_->get<string>("title", ""));   // this is quite a cumbersome way to do this: cleaning needed

  shared_ptr<const Reference> ref = ref_;
  auto m = input->begin();
  for ( ; m != --input->end(); ++m) {
    const std::string title = to_lower((*m)->get<std::string>("title", ""));
    if (title != "molecule") {
      shared_ptr<Method> c = construct_method(title, *m, geom_, ref);
      if (!c) throw runtime_error("unknown method in Hessian");
      c->compute();
      ref = c->conv_to_ref();
    } else {
      geom_ = make_shared<const Geometry>(*geom_, *m);
      if (ref) ref = ref->project_coeff(geom_);
    }
  }
  auto cinput = make_shared<PTree>(**m);
  cinput->put("hessian", true);

  const string method = to_lower(cinput->get<string>("title", ""));

  numhess_ = idata_->get<bool>("numhess", false);
  if (numhess_) {

    cout << "  The numerical Hessian will be computed from gradient calculations" << endl;

    auto displ = std::make_shared<XYZFile>(geom_->natom());
    displ->scale(0.0);

    double plus, minus;
    hess_ = make_shared<Matrix>(3,geom_->natom());

    for (int i = 0; i != geom_->natom(); ++i) {
      for (int j = 0; j != 3; ++j) {
        displ->element(j,i) = 0.0001;
        geom_ = std::make_shared<Geometry>(*geom_, displ);
        geom_->print_atoms();

        auto force_plus = make_shared<Force>(idata_, geom_, ref_);
        force_plus->compute();
        plus = force_plus->grad()->element(j,i);

        displ->element(j,i) = -0.0002;
        geom_ = std::make_shared<Geometry>(*geom_, displ);
        geom_->print_atoms();

        auto force_minus = make_shared<Force>(idata_, geom_, ref_);
        force_minus->compute();
        minus = force_minus->grad()->element(j,i);

        (*hess_)(j,i) = (plus - minus) / 0.0002;      // Hartree / bohr^2

        displ->element(j,i) = 0.0001;
        geom_ = std::make_shared<Geometry>(*geom_, displ);

        displ->element(j,i) = 0.0;
      }
    }

    cout << endl;
    cout << "    * Numerical Hessian matrix";
    for (int j = 0; j != 3; ++j) {
      cout << endl << "      ";
      for (int i = 0; i != geom_->natom(); ++i)
        cout << setw(20) << setprecision(10) << (*hess_)(j,i);
    }
    cout << endl << endl;

  } else {
    cout << "  Analytical Hessian is not currently implemented" << endl;
  }

}
