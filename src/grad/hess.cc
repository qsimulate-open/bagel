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

  numhess_ = idata_->get<bool>("numhess", false);
  numforce_ = idata_->get<bool>("numforce", false);
  if (numhess_)
    if (!numforce_)
      cout << "  The Hessian will be computed with central gradient differences (analytical gradients)" << endl;
    else
      cout << "  The Hessian will be computed with central gradient differences (numerical gradients)" << endl;
  else
    cout << "  Analytical Hessian is not implemented" << endl;

  const string method = to_lower(cinput->get<string>("title", ""));

  if (numhess_) {

  auto displ = std::make_shared<XYZFile>(geom_->natom());
  displ->scale(0.0);
  dx_ = idata_->get<double>("dx", 0.001);
  cout << "  Finite difference size (dx) is " << setprecision(8) << dx_ << " Bohr" << endl;

  muffle_ = make_shared<Muffle>("freq.log");
  hess_ = make_shared<Matrix>(geom_->natom(),3);

  for (int j = 0; j != geom_->natom(); ++j) {  // atom j
    for (int i = 0; i != 3; ++i) { //xyz
      //displace +dx
      displ->element(i,j) = dx_;
      geom_ = std::make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();

      auto plus = make_shared<Force>(idata_, geom_, ref_);
      plus->compute();

      // displace -dx
      displ->element(i,j) = -2.0 * dx_; // undo displacement
      geom_ = std::make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();

      auto minus = make_shared<Force>(idata_, geom_, ref_);
      minus->compute();

      displ->element(i,j) = dx_; //back to original coordinates
      geom_ = make_shared<Geometry>(*geom_, displ);
      displ->element(i,j) = 0.0;

        for (int k = 0; k != geom_->natom(); ++k) {  // atom j
          for (int l = 0; l != 3; ++l) { //xyz

        (*hess_)(k,l) = (plus->grad()->element(k,l) - minus->grad()->element(k,l)) / (2*dx_);
        cout << " hess atom k " << k << " and coord l " << l << "  is " << setprecision(10) << (*hess_)(k,l) << endl;
          }
        }
      }
    }
  muffle_->unmute();

#if 0
  //print hessian
  hess_ = make_shared<Matrix>(3*geom_->natom(),3*geom_->natom());

  cout << endl;
  cout << "    * Numerical Hessian matrix";
    for (int j = 0; j != 3*geom_->natom(); ++j) {
      cout << endl << "      ";
      for (int i = 0; i != 3*geom_->natom(); ++i){
        if (i == j) {
          (*hess_)(i,i) = ((*grad_plus_)(i,i) - (*grad_minus_)(i,i))/(2*dx_);
        } else {
          (*hess_)(i,j) = 0.0;
        }
        cout << setw(20) << setprecision(10) << (*hess_)(i,j);
      }
    }
  cout << endl << endl;
#endif

  }
}
