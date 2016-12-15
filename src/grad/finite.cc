//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: finite.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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

#include <src/grad/gradeval.h>
#include <src/grad/finite.h>
#include <src/util/timer.h>
#include <src/wfn/construct_method.h>

using namespace std;
using namespace bagel;

shared_ptr<GradFile> FiniteGrad::compute() {
  cout << "  Gradient evaluation with respect to " << geom_->natom() * 3 << " DOFs" << endl;
  cout << "  Finite difference size (dx) is " << setprecision(8) << dx_ << " Bohr" << endl;
 
  auto displ = std::make_shared<XYZFile>(geom_->natom());
  auto grad = std::make_shared<GradFile>(geom_->natom());

  double energy_plus, energy_minus;
 
  displ->scale(0.0);
 
  shared_ptr<Method> energy_method;
 
  energy_method = construct_method(method_, idata_, geom_, ref_);
  energy_method->compute();
  ref_ = energy_method->conv_to_ref();
  energy_ = ref_->energy(target_state_);

  // TODO when we start CASSCF calculation with an initial guess from reference geometry,
  // it should be scientifically better than doing the calculation all over again; but it doesn't work
  // Currently the code "works", but in a silly way (for all DOFs, it does SCF again).
  // I didn't delete the part for doing CASSCF calculation starting from reference...
  // this is why the code at current stage is with the total stupidness...
  shared_ptr<const Reference> refgrad_plus;
  shared_ptr<const Reference> refgrad_minus;
 
  cout << "  Reference energy is " << energy_ << endl;
 
  for (int i = 0; i != geom_->natom(); ++i) {
    for (int j = 0; j != 3; ++j) {

      displ->element(j,i) = dx_;
      geom_ = std::make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();

      refgrad_plus = make_shared<Reference>(*ref_, nullptr);
      refgrad_plus = refgrad_plus->project_coeff(geom_);

      energy_method = construct_method(method_, idata_, geom_, refgrad_plus);
      energy_method->compute();
      refgrad_plus = energy_method->conv_to_ref();
      energy_plus = refgrad_plus->energy(target_state_);

      displ->element(j,i) = -2.0 * dx_;
      geom_ = std::make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();
      
      refgrad_minus = make_shared<Reference>(*ref_, nullptr);
      refgrad_minus = refgrad_minus->project_coeff(geom_);

      energy_method = construct_method(method_, idata_, geom_, refgrad_minus);
      energy_method->compute();
      refgrad_minus = energy_method->conv_to_ref();
      energy_minus = refgrad_minus->energy(target_state_);
 
      grad->element(j,i) = (energy_plus - energy_minus) / (dx_ * 2.0);      // Hartree / bohr

      // to the original position

      displ->element(j,i) = dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);

      displ->element(j,i) = 0.0;
    }
  }

  grad->print(": Calculated with finite difference", 0);
  return grad;
}

template<>
shared_ptr<GradFile> FiniteNacm<CASSCF>::compute() {
  cout << "  NACME evaluation with respect to " << geom_->natom() * 3 << " DOFs" << endl;
  cout << "  Finite difference size (dx) is " << setprecision(8) << dx_ << " Bohr(s)" << endl;
  
  auto displ = make_shared<XYZFile>(geom_->natom());
  auto grad_ci = std::make_shared<GradFile>(geom_->natom());
  auto grad = make_shared<GradFile>(geom_->natom());

  displ->scale(0.0);
  
  shared_ptr<Method> energy_method;
  shared_ptr<Dvec> civ_ref = ref_->civectors()->copy();
  int nclosed = ref_->nclosed();
  int nocc = ref_->nocc();
  shared_ptr<const Matrix> acoeff_ref;
  acoeff_ref = make_shared<Matrix>(ref_->coeff()->slice(nclosed, nocc));
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
      displ->element(j,i) = dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();

      refgrad_plus = make_shared<Reference> (*ref_, nullptr);
      refgrad_plus = nullptr;
  
      energy_method = construct_method(method_, idata_, geom_, refgrad_plus);
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

      displ->element(j,i) = -2.0 * dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();

      refgrad_minus = make_shared<Reference> (*ref_, nullptr);
      refgrad_minus = nullptr;
      
      energy_method = construct_method(method_, idata_, geom_, refgrad_minus);
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
      civ_diff->scale(1.0 / (2.0 * dx_));
      acoeff_diff = make_shared<Matrix>(*acoeff_plus - *acoeff_minus);
      acoeff_diff->scale(1.0 / (2.0 * dx_));
  
      displ->element(j,i) = dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);

      displ->element(j,i) = 0.0;

      grad->element(j,i) = civ_ref->data(target_state1_)->dot_product(civ_diff->data(target_state2_));
      grad_ci->element(j,i) = grad->element(j,i);

      auto Uij = make_shared<Matrix>(*acoeff_ref % *Smn * *acoeff_diff);
      for (int ii = 0; ii != norb; ++ii) {
        for (int ij = 0; ij != norb; ++ij) {
          if (ii != ij) {
            for (auto& iter : civ_ref->det()->phia(ii, ij)) {
              size_t iaA = iter.source;
              size_t iaB = iter.target;
              double sign = static_cast<double>(iter.sign);

              for (size_t ib = 0; ib != lenb; ++ib) {                    
                double factor = civ_ref->data(target_state1_)->data(ib+iaB*lenb) * civ_ref->data(target_state2_)->data(ib+iaA*lenb) * sign;
                grad->element(j,i) += factor * (Uij->element(ij, ii) - Uij->element(ii, ij)) * .5;
                if ((i + j * 3) == 0) {
                  gmo->element(ij, ii) += factor * .5;
                  gmo->element(ii, ij) -= factor * .5;
                }

              }
            }
            for (size_t ia = 0; ia != lena; ++ia) {
              for (auto& iter : civ_ref->det()->phib(ii, ij)) {
                size_t ibA = iter.source;
                size_t ibB = iter.target;
                double sign = static_cast<double>(iter.sign);
                double factor = civ_ref->data(target_state1_)->data(ibB+ia*lenb) * civ_ref->data(target_state2_)->data(ibA+ia*lenb) * sign;
                grad->element(j,i) += factor * (Uij->element(ij, ii) - Uij->element(ii, ij)) * .5;
                if ((i + j * 3) == 0) {
                  gmo->element(ij, ii) += factor * .5;
                  gmo->element(ii, ij) -= factor * .5;
                }
              }
            }
          }
        }
      }
    }
  }
  grad_ci->print(": CI term without orbital relaxation, <cJ|d/dXa cI>", 0);
  auto gfin = make_shared<Matrix>(*acoeff_ref * *gmo ^ *acoeff_ref);
  auto grad_basis = make_shared<GradFile>(geom_->natom());
  grad_basis = contract_nacme(nullptr, nullptr, nullptr, nullptr, gfin, /*numerical=*/true);

  grad->print(": CI term, <cJ|d/dXa cI>", 0);
  grad_basis->print(": CSF term, basis set derivative (analytically calculated)", 0);

  *grad += *grad_basis;
  grad->print(": NACME calculated with finite difference", 0);

  return grad;
}
