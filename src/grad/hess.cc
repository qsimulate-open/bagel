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
#include <src/util/atommap.h>
#include <src/util/constants.h>

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
  mw_hess_ = make_shared<Matrix>(3*natom,3*natom);

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
      double energy_ref;

      shared_ptr<Method> energy_method;

      energy_method = construct_method(method, idata_, geom_, ref_);
      energy_method->compute();
      ref_ = energy_method->conv_to_ref();
      energy_ref = ref_->energy(target);

      cout << "  Reference energy is " << setprecision(25)<< energy_ref << endl;

      shared_ptr<const Reference> refgrad_plus;
      shared_ptr<const Reference> refgrad_minus;

      for (int i = 0; i != natom; ++i) {  // atom i
        for (int j = 0; j != 3; ++j) { //xyz
          // displace +dx
          displ->element(j,i) = dx_;
          auto geom_plus = std::make_shared<const Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_plus->print_atoms();

          refgrad_plus = make_shared<Reference> (*ref_, nullptr);
          refgrad_plus = nullptr;

          auto plus = make_shared<FiniteGrad>(method, cinput, geom_plus, refgrad_plus, target, dx_);
          shared_ptr<GradFile> outplus = plus->compute();

          //displace -dx
          displ->element(j,i) = -dx_;
          auto geom_minus = std::make_shared<const Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_minus->print_atoms();

          refgrad_minus = make_shared<Reference> (*ref_, nullptr);
          refgrad_minus = nullptr;

          auto minus = make_shared<FiniteGrad>(method, cinput, geom_minus, refgrad_minus, target, dx_);
          shared_ptr<GradFile> outminus = minus->compute();
          displ->element(j,i) = 0.0;

          for (int k = 0; k != natom; ++k) {  // atom j
            for (int l = 0; l != 3; ++l) { //xyz
              (*hess_)(counter,step) = (outplus->element(l,k) - outminus->element(l,k)) / (2*dx_);
              (*mw_hess_)(counter,step) =  (*hess_)(counter,step) / sqrt(geom_->atoms(i)->averaged_mass() * geom_->atoms(k)->averaged_mass());
              step = step + 1;
            }
          }
          step = 0;
          counter = counter + 1;
        }
      }

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
          displ->element(j,i) = -dx_;
          auto geom_minus = std::make_shared<Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_minus->print_atoms();

          auto minus = make_shared<Force>(idata_, geom_minus, ref_);
          shared_ptr<GradFile> outminus = minus->compute();

          displ->element(j,i) = 0.0;

          for (int k = 0; k != natom; ++k) {  // atom j
            for (int l = 0; l != 3; ++l) { //xyz
              (*hess_)(counter,step) = (outplus->element(l,k) - outminus->element(l,k)) / (2*dx_);
              (*mw_hess_)(counter,step) =  (*hess_)(counter,step) / sqrt(geom_->atoms(i)->averaged_mass() * geom_->atoms(k)->averaged_mass());
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
    hess_->print("Numerical Hessian matrix", 3*natom);

    //print mass weighted hessian
    mw_hess_->print("Mass Weighted Numerical Hessian Matrix", 3*natom);

    //symmetrize mass weighted hessian
    mw_hess_->symmetrize();
    mw_hess_->print("Symmetrized Mass Weighted Numerical Hessian Matrix", 3*natom);
    cout << "    (masses averaged over the natural occurance of isotopes)";
    cout << endl << endl;

    //TODO: Project out rotations and translations

    //diagonalize mass weighted hessian
    // eig(i) in Hartree/bohr^2*amu
    VectorB eig(3*natom);
    mw_hess_->diagonalize(eig);

    // v = sqrt (eig) / (2 pi c )
    for (int i = 0; i != natom; ++i) {
      cout << setw(20) << setprecision(15) << eig(i);
    }
    cout << endl << endl;


    cout << "    * Vibrational Frequencies (wavenumbers)" << endl;
    for (int i = 0; i != 3*natom; ++i) {
      cout << setw(20) << setprecision(7) << sqrt((eig(i) * au2joule__) / amu2kilogram__ ) / (100.0 * au2meter__ * 2.0 * pi__ * csi__);
    }
    cout << endl << endl;
  }
}
