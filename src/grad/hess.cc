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
#include <src/grad/force.h>
#include <src/grad/finite.h>
#include <src/wfn/construct_method.h>
#include <src/grad/gradeval.h>

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
      if (!c) throw runtime_error("unknown method in force");
      c->compute();
      ref = c->conv_to_ref();
    } else {
      geom_ = make_shared<const Geometry>(*geom_, *m);
      if (ref) ref = ref->project_coeff(geom_);
    }
  }
  auto cinput = make_shared<PTree>(**m);
  cinput->put("hessian", true);

  int natom = geom_->natom();
  numhess_ = idata_->get<bool>("numhess", false);
  numforce_ = idata_->get<bool>("numforce", false);
  if (numhess_)
    if (!numforce_)
      cout << "  The Hessian will be computed with central gradient differences (analytical gradients)" << endl;
    else
      cout << "  The Hessian will be computed with central finite difference" << endl;
  else
    cout << "  Analytical Hessian is not implemented" << endl;

  hess_ = make_shared<Matrix>(3*natom,3*natom);

  if (numhess_) {
    dx_ = idata_->get<double>("dx", 0.001);
    cout << "  Finite difference size (dx) is " << setprecision(8) << dx_ << " Bohr" << endl;

    auto displ = std::make_shared<XYZFile>(natom);
    displ->scale(0.0);

    muffle_ = make_shared<Muffle>("freq.log");
    int counter = 0;
    int step = 0;

    if (numforce_) {
      const string method = to_lower(cinput->get<string>("title", ""));
      const int target = idata_->get<int>("target", 0);
      double energy_plus1, energy_plus2, energy_minus1, energy_minus2, energy_ref;

      shared_ptr<Method> energy_method;

      energy_method = construct_method(method, idata_, geom_, ref_);
      energy_method->compute();
      ref_ = energy_method->conv_to_ref();
      energy_ref = ref_->energy(target);

      cout << "  Reference energy is " << energy_ref << endl;

      shared_ptr<const Reference> refgrad_plus1;
      shared_ptr<const Reference> refgrad_plus2;
      shared_ptr<const Reference> refgrad_minus1;
      shared_ptr<const Reference> refgrad_minus2;
      shared_ptr<const Reference> refgrad_plus3;
      shared_ptr<const Reference> refgrad_minus3;

     //compute diagonal elements of hessian
      for (int i = 0; i != natom; ++i) {
        for (int j = 0; j != 3; ++j) {
          // displace +dx
          displ->element(j,i) = dx_;
          auto geom_plus1 = std::make_shared<Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_plus1->print_atoms();

          refgrad_plus1 = make_shared<Reference> (*ref_, nullptr);
          refgrad_plus1 = nullptr;

          energy_method = construct_method(method, idata_, geom_plus1, refgrad_plus1);
          energy_method->compute();
          refgrad_plus1 = energy_method->conv_to_ref();
          energy_plus1 = refgrad_plus1->energy(target);
          displ->element(j,i) = 0;

          // displace +2dx
          displ->element(j,i) = 2.0 * dx_;
          auto geom_plus2 = std::make_shared<Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_plus2->print_atoms();

          refgrad_plus2 = make_shared<Reference> (*ref_, nullptr);
          refgrad_plus2 = nullptr;

          energy_method = construct_method(method, idata_, geom_plus2, refgrad_plus2);
          energy_method->compute();
          refgrad_plus2 = energy_method->conv_to_ref();
          energy_plus2 = refgrad_plus2->energy(target);
          displ->element(j,i) = 0;

          // displace -dx
          displ->element(j,i) = -dx_;
          auto geom_minus1 = std::make_shared<Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_minus1->print_atoms();

          refgrad_minus1 = make_shared<Reference> (*ref_, nullptr);
          refgrad_minus1 = nullptr;

          energy_method = construct_method(method, idata_, geom_minus1, refgrad_minus1);
          energy_method->compute();
          refgrad_minus1 = energy_method->conv_to_ref();
          energy_minus1 = refgrad_minus1->energy(target);
          displ->element(j,i) = 0;

          // displace -2dx
          displ->element(j,i) = -2.0 * dx_;
          auto geom_minus2 = std::make_shared<Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_minus2->print_atoms();

          refgrad_minus2 = make_shared<Reference> (*ref_, nullptr);
          refgrad_minus2 = nullptr;

          energy_method = construct_method(method, idata_, geom_minus2, refgrad_minus2);
          energy_method->compute();
          refgrad_minus2 = energy_method->conv_to_ref();
          energy_minus2 = refgrad_minus2->energy(target);
          displ->element(j,i) = 0;

          (*hess_)(counter,counter) = (-energy_plus2 + (16 * energy_plus1) - (30 * energy_ref) + (16 * energy_minus1) - energy_minus2) / ( 12 * dx_ * dx_);
          counter = counter + 1;
        }
      }

      counter = 0;
      step = 0;

#if 0
      for (int i = 0; i != natom; ++i) {  // atom i
        for (int j = 0; j != 3; ++j) { //xyz

          //displace +dx
          displ->element(j,i) = dx_;
          auto geom_plus = std::make_shared<const Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_plus->print_atoms();

          refgrad_plus3 = make_shared<Reference> (*ref_, nullptr);
          refgrad_plus3 = nullptr;

          auto plus = make_shared<FiniteGrad>(method, cinput, geom_plus, refgrad_plus3, target, dx_);
          shared_ptr<GradFile> outplus = plus->compute();

          // displace -dx
          displ->element(j,i) = -dx_; // undo displacement
          auto geom_minus = std::make_shared<Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_minus->print_atoms();

          refgrad_minus3 = make_shared<Reference> (*ref_, nullptr);
          refgrad_minus3 = nullptr;

          auto minus = make_shared<Force>(method, cinput, geom_minus, refgrad_minus3, target, dx_);
          shared_ptr<GradFile> outminus = minus->compute();

          displ->element(j,i) = 0.0;
          if (counter != step) {
            for (int k = 0; k != natom; ++k) {  // atom j
              for (int l = 0; l != 3; ++l) { //xyz
                //storing hessian in a really silly way. Need to improve
                (*hess_)(counter,step) = (outplus->element(l,k) - outminus->element(l,k)) / (2*dx_);
                step = step + 1;
              }
            }
          }
          step = 0;
          counter = counter + 1;
        }
      }
#endif



    } else {  //finite difference with analytical gradients

      for (int i = 0; i != natom; ++i) {  // atom i
        for (int j = 0; j != 3; ++j) { //xyz

          //displace +dx
          displ->element(j,i) = dx_;
          auto geom_plus = std::make_shared<const Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_plus->print_atoms();

          auto plus = make_shared<Force>(idata_, geom_plus, ref_);
          shared_ptr<GradFile> outplus = plus->compute();

          // displace -dx
          displ->element(j,i) = -dx_; // undo displacement
          auto geom_minus = std::make_shared<Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_minus->print_atoms();

          auto minus = make_shared<Force>(idata_, geom_minus, ref_);
          shared_ptr<GradFile> outminus = minus->compute();

          displ->element(j,i) = 0.0;

          for (int k = 0; k != natom; ++k) {  // atom j
            for (int l = 0; l != 3; ++l) { //xyz
              //storing hessian in a really silly way. Need to improve
              (*hess_)(counter,step) = (outplus->element(l,k) - outminus->element(l,k)) / (2*dx_);
              step = step + 1;
            }
          }
          step = 0;
          counter = counter + 1;
        }
      }
    }
    muffle_->unmute();

    //print hessian
    cout << endl;
    cout << "    * Numerical Hessian matrix";
    for (int i = 0; i != 3*natom; ++i) {
      cout << endl << "      ";
      for (int j = 0; j != 3*natom; ++j){
        cout << setw(20) << setprecision(10) << (*hess_)(j,i);
      }
    }
    cout << endl << endl;

    //print symmetrized hessian
#if 1
    cout << endl;
    cout << "    * Symmetrized Numerical Hessian matrix";
    for (int j = 0; j != 3*geom_->natom(); ++j) {
      cout << endl << "      ";
      for (int i = 0; i != 3*geom_->natom(); ++i){
        cout << setw(20) << setprecision(10) << ((*hess_)(i,j) + (*hess_)(j,i)) / 2.0 ;
      }
    }
    cout << endl << endl;
#endif

  }
}
