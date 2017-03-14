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
#include <src/util/timer.h>
#include <src/prop/multipole.h>

using namespace std;
using namespace bagel;

Hess::Hess(shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : idata_(idata), geom_(g), ref_(r) {

}
void Hess::compute() {
  auto input = idata_->get_child("method");
  const string jobtitle = to_lower(idata_->get<string>("title", ""));   // this is quite a cumbersome way to do this: cleaning needed
  Timer timer;

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
  auto cartesian = make_shared<Matrix>(3,3*natom); //matrix of dmu/dR

  if (numhess_) {
    dx_ = idata_->get<double>("dx", 0.001);
    cout << "  Finite difference size (dx) is " << setprecision(8) << dx_ << " Bohr" << endl;

    std::vector<double> dipole_minus;
    std::vector<double> dipole_plus;

    auto displ = std::make_shared<XYZFile>(natom);
    displ->scale(0.0);

    muffle_ = make_shared<Muffle>("freq.log");
    int counter = 0;
    int step = 0;
    if (numforce_) {
      muffle_->mute();
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

          displ->element(j,i) = dx_;
          auto geom_plus = std::make_shared<const Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_plus->print_atoms();

          refgrad_plus = make_shared<Reference> (*ref_, nullptr);
          refgrad_plus = nullptr;

          auto plus = make_shared<FiniteGrad>(method, cinput, geom_plus, refgrad_plus, target, dx_);
          shared_ptr<GradFile> outplus = plus->compute();

//          dipole_plus = plus->ref()->dipole();

          //displace -dx
          displ->element(j,i) = -dx_;
          auto geom_minus = std::make_shared<const Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_minus->print_atoms();

          refgrad_minus = make_shared<Reference> (*ref_, nullptr);
          refgrad_minus = nullptr;

          auto minus = make_shared<FiniteGrad>(method, cinput, geom_minus, refgrad_minus, target, dx_);
          shared_ptr<GradFile> outminus = minus->compute();

//          dipole_minus = minus->ref()->dipole();

          displ->element(j,i) = 0.0;

          for (int k = 0; k != natom; ++k) {  // atom j
            for (int l = 0; l != 3; ++l) { //xyz
              (*hess_)(counter,step) = (outplus->element(l,k) - outminus->element(l,k)) / (2*dx_);
              (*mw_hess_)(counter,step) =  (*hess_)(counter,step) / sqrt(geom_->atoms(i)->averaged_mass() * geom_->atoms(k)->averaged_mass());
//              (*cartesian)(l,counter) = (dipole_plus[l] - dipole_minus[l]) / (2*dx_);
              step = step + 1;
            }
          }
          step = 0;
          counter = counter + 1;
          muffle_->unmute();
        }
      }

    } else {  //finite difference with analytical gradients
      for (int i = 0; i != natom; ++i) {  // atom i
        for (int j = 0; j != 3; ++j) { //xyz
          // displace +dx
          muffle_->mute();

          //displace +dx
          displ->element(j,i) = dx_;
          auto geom_plus = std::make_shared<const Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_plus->print_atoms();
          if (ref_)
            ref_ = ref_->project_coeff(geom_plus);

          auto plus = make_shared<Force>(idata_, geom_plus, ref_);
          shared_ptr<GradFile> outplus = plus->compute();
          dipole_plus = plus->force_dipole();

          // displace -dx
          displ->element(j,i) = -dx_;
          auto geom_minus = std::make_shared<Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);
          geom_minus->print_atoms();
          if (ref_)
            ref_ = ref_->project_coeff(geom_minus);

          auto minus = make_shared<Force>(idata_, geom_minus, ref_);
          shared_ptr<GradFile> outminus = minus->compute();
          dipole_minus = minus->force_dipole();

          displ->element(j,i) = 0.0;

          for (int k = 0; k != natom; ++k) {  // atom j
            for (int l = 0; l != 3; ++l) { //xyz
              (*hess_)(counter,step) = (outplus->element(l,k) - outminus->element(l,k)) / (2*dx_);
              (*mw_hess_)(counter,step) =  (*hess_)(counter,step) / sqrt(geom_->atoms(i)->averaged_mass() * geom_->atoms(k)->averaged_mass());
              (*cartesian)(l,counter) = (dipole_plus[l] - dipole_minus[l]) / (2*dx_);
              step = step + 1;
            }
          }
          step = 0;
          counter = counter + 1;
          muffle_->unmute();
          stringstream ss; ss << "Hessian evaluation (" << setw(2) << i*3+j+1 << " / " << natom * 3 << ")";
          timer.tick_print(ss.str());
        }
      }
    }

    //symmetrize mass weighted hessian
    hess_->print("Hessian");
    mw_hess_->symmetrize();
    cout << "    (masses averaged over the natural occurance of isotopes)";
    cout << endl << endl;
    mw_hess_->print("Symmetrized Mass Weighted Hessian", 3*natom);

    // calculate center of mass
    VectorB num(3); // values needed to calc center of mass. mi*xi, mi*yi, mi*zi, and total mass
    VectorB center(3);
    double total_mass = 0;
    // compute center of mass
    for (int i = 0; i!= natom; ++i) {
      for (int j = 0; j != 3; ++j) {
        num(j) += geom_->atoms(i)->averaged_mass() * geom_->atoms(i)->position(j);
      }
      total_mass += geom_->atoms(i)->averaged_mass();
    }

    for (int j = 0; j!= 3; ++j) {
      center(j) = num(j)/total_mass;
    }

    // shift center coordinate system to center of mass
    for (int i = 0; i != natom; ++i) {
      for (int j = 0; j != 3; ++j) {
        displ->element(j,i) = -1.0 * center(j);
      }
    }
    auto geom_center = std::make_shared<Geometry>(*geom_, displ, make_shared<const PTree>(), false, false);

    cout << "    * Projecting out Translational Degrees of Freedom " << endl;
    auto identity =  make_shared<Matrix>(3*natom,3*natom);
    for (int ist = 0; ist!= 3*natom; ++ist) {
      for (int jst = 0; jst != 3*natom; ++jst) {
        if (jst == ist)
          (*identity)(jst,ist) = 1;
      }
    }
    // 3N by 3 matrix P_trans
    auto proj_trans =  make_shared<Matrix>(3,3*natom);
    for (int i = 0; i!= natom; ++i) {
      for (int j = 0; j != 3; ++ j) {
        (*proj_trans)(j, 3*i+j) = sqrt(geom_->atoms(i)->averaged_mass());
      }
    }

    cout << "    * Projecting out Rotational Degrees of Freedom " << endl;

    // 3N by 3 matrix P of orthogonal vectors about the XYZ axes
    auto proj_rot = make_shared<Matrix>(3,3*natom);
    for (int i = 0; i != natom; ++i) {
      (*proj_rot)(0,3*i) = 0;
      (*proj_rot)(0,3*i+1) = -1.0 * sqrt(geom_->atoms(i)->averaged_mass()) * (geom_center->atoms(i)->position(2));
      (*proj_rot)(0,3*i+2) = sqrt(geom_->atoms(i)->averaged_mass()) * geom_center->atoms(i)->position(1);

      (*proj_rot)(1,3*i) =  sqrt(geom_->atoms(i)->averaged_mass()) * geom_center->atoms(i)->position(2);
      (*proj_rot)(1,3*i+1) = 0;
      (*proj_rot)(1,3*i+2) = -1.0 * sqrt(geom_->atoms(i)->averaged_mass()) * (geom_center->atoms(i)->position(0));

      (*proj_rot)(2,3*i) = -1.0 * sqrt(geom_->atoms(i)->averaged_mass()) * (geom_center->atoms(i)->position(1));
      (*proj_rot)(2,3*i+1) = sqrt(geom_->atoms(i)->averaged_mass()) * geom_center->atoms(i)->position(0);
      (*proj_rot)(2,3*i+2) = 0;
    }

    auto proj_total = make_shared<Matrix>(6, 3*natom);
    for (int i = 0; i != 3*natom; ++i) {
      (*proj_total)(0, i) = (*proj_trans)(0,i);
      (*proj_total)(1, i) = (*proj_trans)(1,i);
      (*proj_total)(2, i) = (*proj_trans)(2,i);
      (*proj_total)(3, i) = (*proj_rot)(0,i);
      (*proj_total)(4, i) = (*proj_rot)(1,i);
      (*proj_total)(5, i) = (*proj_rot)(2,i);
    }

    //normalize the set of six orthogonal vectors
    double trans1 = 0;
    double trans2 = 0;
    double trans3 = 0;
    double rot1 = 0;
    double rot2 = 0;
    double rot3 = 0;

    // sum of the square of each vector
    for (int i = 0; i != 3*natom; ++i) {
      trans1 += (*proj_total)(0,i) * (*proj_total)(0,i);
      trans2 += (*proj_total)(1,i) * (*proj_total)(1,i);
      trans3 += (*proj_total)(2,i) * (*proj_total)(2,i);
      rot1 += (*proj_total)(3,i) * (*proj_total)(3,i);
      rot2 += (*proj_total)(4,i) * (*proj_total)(4,i);
      rot3 += (*proj_total)(5,i) * (*proj_total)(5,i);
    }

    auto proj_norm = make_shared<Matrix>(6, 3*natom);
    for (int i = 0; i != 3*natom; ++i) {
      (*proj_norm)(0,i) = (*proj_total)(0,i) / sqrt(trans1);
      (*proj_norm)(1,i) = (*proj_total)(1,i) / sqrt(trans2);
      (*proj_norm)(2,i) = (*proj_total)(2,i) / sqrt(trans3);
      (*proj_norm)(3,i) = (*proj_total)(3,i) / sqrt(rot1);
      (*proj_norm)(4,i) = (*proj_total)(4,i) / sqrt(rot2);
      if (rot3 != 0)
        (*proj_norm)(5,i) = (*proj_total)(5,i) / sqrt(rot3);
    }

    auto proj_hess_ = make_shared<Matrix>(3*natom,3*natom);
    *proj_hess_ = (*identity - *proj_norm % *proj_norm) % *mw_hess_ * (*identity - *proj_norm % *proj_norm);

    //diagonalize hessian
    // eig(i) in Hartree/bohr^2*amu
    VectorB eig(3*natom);
    proj_hess_->diagonalize(eig);

    cout << endl << " Mass Weighted Hessian Eigenvalues" << endl; // units Hartree/bohr^2*amu
    for (int i = 0; i != 3*natom; ++i ) {
      cout << setw(10) << setprecision(5) << eig(i);
    }
    cout << endl;

    proj_hess_->print("Mass Weighted Hessian Eigenvectors", 3*natom);

    //convert mw eigenvectors to normalized cartesian modes
    auto eigvec_cart = make_shared<Matrix>(3*natom,3*natom);

    counter = 0.0;
    step = 0.0;
    for (int i = 0; i != natom; ++i) {
      for (int j = 0; j != 3; ++j) {
        for (int k = 0; k != natom; ++k) {
          for (int l = 0; l != 3; ++l) {
            (*eigvec_cart)(step,counter) =  proj_hess_->element(step,counter) / sqrt(geom_->atoms(k)->averaged_mass());
            step = step + 1;
          }
        }
        step = 0;
        counter = counter + 1;
      }
    }

    // calculate IR intensities: (N*pi/3c^2) * (dmu/dq)^2
    auto normal = make_shared<Matrix>(3,3*natom); //matrix of dmu/dQ
    *normal =*cartesian * *eigvec_cart;  // units (e bohr / bohr) * 1/sqrt(amu)

    VectorB dmudq2(3*natom); //square of the dipole derivatives in hartree * bohr / amu
    for (int i = 0; i != 3*natom; ++i) { // total for each mode:
        dmudq2(i) =((*normal)(0,i)*(*normal)(0,i) + (*normal)(1,i)*(*normal)(1,i) + (*normal)(2,i)*(*normal)(2,i));
    }

    cout << "    * Vibrational frequencies, IR intensities, and corresponding cartesian eigenvectors" << endl << endl;

    int len_n = eigvec_cart->ndim();
    int len_m = eigvec_cart->mdim();

    for (int i = 0; i != len_m/6; ++i) {
      cout << setw(15) << " ";
      for (int k = 0; k != 6; ++k)
        cout << setw(20) << i * 6 + k;
      cout << endl << setw(15) << "Freq (cm**-1)";
      for (int k = 0; k != 6; ++k)
        cout << setw(20) << setprecision(2) << (fabs(eig(i*6+k)) > 1.0e-6 ? (eig(i*6+k) > 0.0 ? sqrt((eig(i*6+k) * au2joule__) / amu2kilogram__) / (100.0 * au2meter__ * 2.0 * pi__ * csi__) : -sqrt((-eig(i*6+k) * au2joule__) / amu2kilogram__) / (100.0 * au2meter__ * 2.0 * pi__ * csi__)) : 0) ;
      cout << endl;

      cout << endl << setw(15) << "IR Int. (km/mol)";
      for (int k = 0; k != 6; ++k)
        cout << setw(20) << setprecision(2) <<(fabs(eig(i*6+k)) > 1.0e-6 ? (eig(i*6+k) > 0.0 ? (((avogadro__ * pi__ *au2meter__ * au2joule__)/ (3 * 1000 * csi__ * csi__ * amu2kilogram__)) * dmudq2(i*6+k)) : (((avogadro__ * pi__ *au2meter__ * au2joule__)/ (3 * 1000 * csi__ * csi__ * amu2kilogram__)) * dmudq2(i*6+k))) : 0) ;
      cout << endl << endl;

      for (int j = 0; j != len_n; ++j) {
        cout << setw(15) << j;
        for (int k = 0; k != 6; ++k)
          cout << setw(20) << setprecision(5) << eigvec_cart->element(j, i*6+k);
        cout << endl;
      }
      cout << endl;
    }

    if (len_m%6) {
      int i = len_m/6;
      cout << setw(15) << " ";
      for (int k = 0; k != len_m%6; ++k)
        cout << setw(20) << i * 6 + k;
      cout << endl << setw(15) << "Freq (cm**-1)";
      for (int k = 0; k != len_m%6; ++k)
        cout << setw(20) << setprecision(2) << (fabs(eig(i*6+k)) > 1.0e-6 ? (eig(i*6+k) > 0.0 ? sqrt((eig(i*6+k) * au2joule__) / amu2kilogram__) / (100.0 * au2meter__ * 2.0 * pi__ * csi__) : -sqrt((-eig(i*6+k) * au2joule__) / amu2kilogram__) / (100.0 * au2meter__ * 2.0 * pi__ * csi__)) : 0) ;
      cout << endl ;

      cout << endl << setw(15) << "IR Int. (a.u.)";
      for (int k = 0; k != len_m%6; ++k)
        cout << setw(20) << setprecision(2) <<(fabs(eig(i*6+k)) > 1.0e-6 ? (eig(i*6+k) > 0.0 ? (((avogadro__ * pi__ *au2meter__ * au2joule__)/ (3 * 1000 * csi__ * csi__ * amu2kilogram__)) * dmudq2(i*6+k)) : (((avogadro__ * pi__ *au2meter__ * au2joule__)/ (3 * 1000 * csi__ * csi__ * amu2kilogram__)) * dmudq2(i*6+k))) : 0) ;
      cout << endl << endl;

      for (int j = 0; j != len_n; ++j) {
        cout << setw(15) << j;
        for (int k = 0; k != len_m%6; ++k)
          cout << setw(20) << setprecision(5) << eigvec_cart->element(j, i*6+k);
        cout << endl;
      }
      cout << endl;
    }

  }
}
