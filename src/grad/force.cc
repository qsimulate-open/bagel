//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: force.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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
#include <src/grad/force.h>
#include <src/grad/gradeval.h>
#include <src/wfn/construct_method.h>

using namespace std;
using namespace bagel;

Force::Force(shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : idata_(idata), geom_(g), ref_(r) {

}


void Force::compute() {
  auto input = idata_->get_child("method");
  const int target = idata_->get<int>("target", 0);
  const string jobtitle = to_lower(idata_->get<string>("title", ""));   // this is quite a cumbersome way to do this: cleaning needed
  const int target2= idata_->get<int>("target2", 1);

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
  cinput->put("gradient", true);

  numerical_ = idata_->get<bool>("numerical", false);
  if (numerical_ == true)
    cout << "  The gradients will be computed with finite difference" << endl;
  else
    cout << "  The gradients will be computed analytically" << endl;


  const string method = to_lower(cinput->get<string>("title", ""));

  if (numerical_ == false) {
    if (method == "uhf") {
 
      auto force = make_shared<GradEval<UHF>>(cinput, geom_, ref_, target);
      force->compute();
 
    } else if (method == "rohf") {
 
      auto force = make_shared<GradEval<ROHF>>(cinput, geom_, ref_, target);
      force->compute();
 
    } else if (method == "hf") {
 
      auto force = make_shared<GradEval<RHF>>(cinput, geom_, ref_, target);
      force->compute();
 
    } else if (method == "ks") {
 
      auto force = make_shared<GradEval<KS>>(cinput, geom_, ref_, target);
      force->compute();
 
    } else if (method == "dhf") {
 
      auto force = make_shared<GradEval<Dirac>>(cinput, geom_, ref_, target);
      force->compute();
 
    } else if (method == "mp2") {
 
      auto force = make_shared<GradEval<MP2Grad>>(cinput, geom_, ref_, target);
      force->compute();
 
    } else if (method == "casscf" && jobtitle == "nacme") {
      
      auto force = make_shared<NacmEval<CASSCF>>(cinput, geom_, ref_, target, target2);
      force->compute();
 
    } else if (method == "casscf") {
 
      auto force = make_shared<GradEval<CASSCF>>(cinput, geom_, ref_, target);
      force->compute();
 
    } else if (method == "caspt2") {
 
      auto force = make_shared<GradEval<CASPT2Grad>>(cinput, geom_, ref_, target);
      force->compute();
 
    } else {
 
        numerical_ = true;
        cout << "  It seems like no analytical gradient method available; moving to finite difference " << endl;
 
    }
  }

  if (numerical_ == true) {

    const int diffsize = idata_->get<int>("diffsize", 3);
    const double dx = pow(0.1, static_cast<double>(diffsize));

    if (jobtitle == "force") {

      // Refer to the optimization code, opt/opt.h 
      cout << "  Gradient evaluation with respect to " << geom_->natom() * 3 << " DOFs" << endl;
      cout << "  Finite difference size (dx) is " << setprecision(8) << dx << " Bohr" << endl;
     
      auto displ = std::make_shared<XYZFile>(geom_->natom());
      auto gradient = std::make_shared<GradFile>(geom_->natom());

      double energy_plus, energy_minus, energy_ref;
     
      displ->scale(0.0);
     
      shared_ptr<Method> energy_method;
     
      energy_method = construct_method(method, cinput, geom_, ref);
      energy_method->compute();
      ref = energy_method->conv_to_ref();
      energy_ref = ref->energy(target);

      // TODO when we start CASSCF calculation with an initial guess from reference geometry,
      // it should be scientifically better than doing the calculation all over again; but it doesn't work
      // Currently the code "works", but in a silly way (for all DOFs, it does SCF again).
      // I didn't delete the part for doing CASSCF calculation starting from reference...
      // this is why the code at current stage is with the total stupidness...
      shared_ptr<const Reference> refgrad_plus;
      shared_ptr<const Reference> refgrad_minus;
     
      cout << "  Reference energy is " << energy_ref << endl;
     
      for (int i = 0; i != geom_->natom(); ++i) {
        for (int j = 0; j != 3; ++j) {

          displ->element(j,i) = dx;
          geom_ = std::make_shared<Geometry>(*geom_, displ);
          geom_->print_atoms();

          refgrad_plus = make_shared<Reference> (*ref, nullptr);
          refgrad_plus = nullptr;

          energy_method = construct_method(method, cinput, geom_, refgrad_plus);
          energy_method->compute();
          refgrad_plus = energy_method->conv_to_ref();
          energy_plus = refgrad_plus->energy(target);

          displ->element(j,i) = -2.0 * dx;
          geom_ = std::make_shared<Geometry>(*geom_, displ);
          geom_->print_atoms();
          
          refgrad_minus = make_shared<Reference> (*ref, nullptr);
          refgrad_minus = nullptr;

          energy_method = construct_method(method, cinput, geom_, refgrad_minus);
          energy_method->compute();
          refgrad_minus = energy_method->conv_to_ref();
          energy_minus = refgrad_minus->energy(target);
     
          gradient->element(j,i) = (energy_plus - energy_minus) / (dx * 2.0);      // Hartree / bohr

          // to the original position

          displ->element(j,i) = dx;
          geom_ = std::make_shared<Geometry>(*geom_, displ);

          displ->element(j,i) = 0.0;
        }
      }

      gradient->print(": Calculated with finite difference", 0);

    } else if (jobtitle == "nacme") {

      if (method != "casscf") 
        throw logic_error ("  Numerical NACME with the methods != CASSCF not implemented yet");

      cout << "  NACME evaluation with respect to " << geom_->natom() * 3 << " DOFs" << endl;
      cout << "  Finite difference size (dx) is " << setprecision(8) << dx << " Bohr(s)" << endl;
     
      auto displ = std::make_shared<XYZFile>(geom_->natom());
      auto gradient_ci = std::make_shared<GradFile>(geom_->natom());
      auto gradient = std::make_shared<GradFile>(geom_->natom());

      displ->scale(0.0);
      
      auto force = make_shared<NacmEval<CASSCF>>(cinput, geom_, ref, target, target2);
      ref = nullptr;
     
      shared_ptr<Method> energy_method;
      energy_method = construct_method(method, cinput, geom_, ref);
      energy_method->compute();
      ref = energy_method->conv_to_ref();
      shared_ptr<const Dvec> civ_ref = ref->civectors();
      int nclosed = ref->nclosed();
      int nocc = ref->nocc();
      shared_ptr<const Matrix> acoeff_ref;
      acoeff_ref = make_shared<Matrix>(ref->coeff()->slice(nclosed, nocc));
      acoeff_ref->print("acoeff_ref = ");
      civ_ref->print (/*sort=*/false);
      
      shared_ptr<const Reference> refgrad_plus;
      shared_ptr<const Reference> refgrad_minus;
      shared_ptr<Dvec> civ_plus;
      shared_ptr<Dvec> civ_minus;
      shared_ptr<Dvec> civ_diff;
      shared_ptr<Matrix> acoeff_plus;
      shared_ptr<Matrix> acoeff_minus;
      shared_ptr<Matrix> acoeff_diff;

      auto Smn = make_shared<Overlap>(geom_);
      const int norb = civ_ref->det()->norb();
      const int lena = civ_ref->det()->lena();
      const int lenb = civ_ref->det()->lenb();
      auto gmo = make_shared<Matrix>(norb, norb);
      gmo->zero();
      
      assert(norb==(nocc-nclosed));
     
      for (int i = 0; i != geom_->natom(); ++i) {
        for (int j = 0; j != 3; ++j) {
          displ->element(j,i) = dx;
          geom_ = make_shared<Geometry>(*geom_, displ);
          geom_->print_atoms();

          refgrad_plus = make_shared<Reference> (*ref, nullptr);
          refgrad_plus = nullptr;
     
          energy_method = construct_method(method, cinput, geom_, refgrad_plus);
          energy_method->compute();
          refgrad_plus = energy_method->conv_to_ref();
          acoeff_plus = make_shared<Matrix>(refgrad_plus->coeff()->slice(nclosed, nocc));
          for (int im = 0; im != acoeff_ref->mdim(); ++im) {
            double dmatch = blas::dot_product(acoeff_ref->element_ptr(0, im), acoeff_ref->ndim(), acoeff_plus->element_ptr(0,im));
            if (dmatch < 0.0) 
              blas::scale_n(-1.0, acoeff_plus->element_ptr(0, im), acoeff_ref->ndim());
          }

          civ_plus = refgrad_plus->civectors()->copy();
          civ_plus->match(civ_ref);

          displ->element(j,i) = -2.0 * dx;
          geom_ = make_shared<Geometry>(*geom_, displ);
          geom_->print_atoms();

          refgrad_minus = make_shared<Reference> (*ref, nullptr);
          refgrad_minus = nullptr;
          
          energy_method = construct_method(method, cinput, geom_, refgrad_minus);
          energy_method->compute();
          refgrad_minus = energy_method->conv_to_ref();
          acoeff_minus = make_shared<Matrix>(refgrad_minus->coeff()->slice(nclosed, nocc));
          for (int im = 0; im != acoeff_ref->mdim(); ++im) {
            double dmatch = blas::dot_product(acoeff_ref->element_ptr(0, im), acoeff_ref->ndim(), acoeff_minus->element_ptr(0,im));
            if (dmatch < 0.0) 
              blas::scale_n(-1.0, acoeff_minus->element_ptr(0, im), acoeff_ref->ndim());
          }

          civ_minus = refgrad_minus->civectors()->copy();
          civ_minus->match(civ_ref);

          civ_diff = civ_plus->copy();
          *civ_diff -= *civ_minus;
          civ_diff->scale(1.0 / (2.0 * dx));
          civ_diff->print(/*sort=*/false);
          acoeff_diff = make_shared<Matrix>(*acoeff_plus - *acoeff_minus);
          acoeff_diff->scale(1.0 / (2.0 * dx));
     
          displ->element(j,i) = dx;
          geom_ = std::make_shared<Geometry>(*geom_, displ);

          displ->element(j,i) = 0.0;

          gradient->element(j,i) = civ_ref->data(target)->dot_product(civ_diff->data(target2));
          gradient_ci->element(j,i) = gradient->element(j,i);

          auto Kfactor = make_shared<Matrix>(*acoeff_ref % *Smn * *acoeff_diff);
          for (int ii = 0; ii != norb; ++ii) {
            for (int ij = 0; ij != norb; ++ij) {
              if (ii != ij) {
                for (auto& iter : civ_ref->det()->phia(ii, ij)) {
                  size_t iaA = iter.source;
                  size_t iaB = iter.target;
                  double sign = static_cast<double>(iter.sign);

                  for (size_t ib = 0; ib != lenb; ++ib) {                    
                    double factor = civ_ref->data(target)->data(ib+iaB*lenb) * civ_ref->data(target2)->data(ib+iaA*lenb) * sign;
                    gradient->element(j,i) += factor * Kfactor->element(ij, ii);
                    if ((i + j * 3) == 0)
                      gmo->element(ij, ii) += factor;
                  }
                }
                for (size_t ia = 0; ia != lena; ++ia) {
                  for (auto& iter : civ_ref->det()->phib(ii, ij)) {
                    size_t ibA = iter.source;
                    size_t ibB = iter.target;
                    double sign = static_cast<double>(iter.sign);
                    double factor = civ_ref->data(target)->data(ibB+ia*lenb) * civ_ref->data(target2)->data(ibA+ia*lenb) * sign;
                    gradient->element(j,i) += factor * Kfactor->element(ij, ii);
                    if ((i + j * 3) == 0)
                      gmo->element(ij, ii) += factor;
                  }
                }
              }
            }
          }

        }
      }
      auto gfin = make_shared<Matrix>(*acoeff_ref * *gmo ^ *acoeff_ref);
      auto grad_basis = std::make_shared<GradFile>(geom_->natom());
      grad_basis = force->contract_nacme(nullptr, nullptr, nullptr, nullptr, gfin, /*numerical=*/true);

      gradient_ci->print(": CI term", 0);
      gradient->print(": All symmetric terms", 0);
      grad_basis->print(": Basis set derivative (analytically calculated)", 0);

      *gradient += *grad_basis;
      gradient->print(": NACME calculated with finite difference", 0);

    }
  }
}
