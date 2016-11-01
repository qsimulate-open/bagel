//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: caspt2grad.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <bagel_config.h>

#include <src/scf/hf/fock.h>
#include <src/grad/cpcasscf.h>
#include <src/grad/gradeval.h>
#include <src/multi/casscf/cassecond.h>
#include <src/multi/casscf/casnoopt.h>
#include <src/multi/casscf/qvec.h>
#include <src/smith/smith.h>
#include <src/smith/caspt2grad.h>
#include <src/prop/multipole.h>
#include <src/prop/hyperfine.h>


using namespace std;
using namespace bagel;

CASPT2Ener::CASPT2Ener(shared_ptr<const PTree> inp, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
  : Method(inp, geom, ref) {
#ifdef COMPILE_SMITH
  Timer timer;

  // compute CASSCF first
  if (inp->get<string>("algorithm", "") != "noopt") {
    auto cas = make_shared<CASSecond>(inp, geom, ref);
    cas->compute();
    ref_ = cas->conv_to_ref();
    fci_ = cas->fci();
    thresh_ = cas->thresh();
  } else {
    auto cas = make_shared<CASNoopt>(inp, geom, ref);
    cas->compute();
    ref_ = cas->conv_to_ref();
    fci_ = cas->fci();
    thresh_ = cas->thresh();
  }

  // gradient/property calculation
  do_hyperfine_ = inp->get<bool>("hyperfine", false);

  timer.tick_print("Reference calculation");

  cout << endl << "  === DF-CASPT2 calculation ===" << endl << endl;
#else
  throw logic_error("CASPT2 gradients require SMITH-generated code. Please compile BAGEL with --enable-smith");
#endif
}


// compute smith and set rdms and ci deriv to a member
void CASPT2Ener::compute() {
#ifdef COMPILE_SMITH
//  const int nclosed = ref_->nclosed();
//  const int nact = ref_->nact();
//  const int nocc = ref_->nocc();

  // construct SMITH here
  shared_ptr<PTree> smithinput = idata_->get_child("smith");
  smithinput->put("_grad", false);
  smithinput->put("_hyperfine", do_hyperfine_);
  auto smith = make_shared<Smith>(smithinput, ref_->geom(), ref_);
  smith->compute();
  
  msrot_   = smith->msrot();
  energy_  = smith->algo()->energyvec();
  coeff_   = smith->coeff();
#endif
}

template<>
shared_ptr<GradFile> FiniteNacm<CASPT2Ener>::compute() {
#ifdef COMPILE_SMITH
  cout << "  NACME evaluation with respect to " << geom_->natom() * 3 << " DOFs" << endl;
  cout << "  Finite difference size (dx) is " << setprecision(8) << dx_ << " Bohr(s)" << endl;
  
  auto displ = make_shared<XYZFile>(geom_->natom());
  auto grad_ci = std::make_shared<GradFile>(geom_->natom());
  auto grad = make_shared<GradFile>(geom_->natom());

  displ->scale(0.0);
  
  shared_ptr<Dvec> civ_ref = ref_->civectors()->copy();
  int nclosed = ref_->nclosed();
  int nocc = ref_->nocc();
  shared_ptr<const Matrix> acoeff_ref;
  acoeff_ref = make_shared<Matrix>(task_->coeff()->slice(nclosed, nocc));
  acoeff_ref->print("acoeff_ref = ");

  civ_ref->rotate (task_->msrot());
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

  auto idata_out = std::make_shared<PTree>(*idata_);
  idata_out->put("_target", target_state1_);
  idata_out->put("_target2", target_state2_);
  
  for (int i = 0; i != geom_->natom(); ++i) {
    for (int j = 0; j != 3; ++j) {
      displ->element(j,i) = dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();

      refgrad_plus = make_shared<Reference> (*ref_, nullptr);
      refgrad_plus = nullptr;

      task_ = std::make_shared<CASPT2Ener>(idata_out, geom_, refgrad_plus);
      task_->compute();
      refgrad_plus  = task_->conv_to_ref();

      acoeff_plus = make_shared<Matrix>(task_->coeff()->slice(nclosed, nocc));
      for (int im = 0; im != acoeff_ref->mdim(); ++im) {
        double dmatch = blas::dot_product(acoeff_ref->element_ptr(0, im), acoeff_ref->ndim(), acoeff_plus->element_ptr(0,im));
        if (dmatch < 0.0) 
          blas::scale_n(-1.0, acoeff_plus->element_ptr(0, im), acoeff_ref->ndim());
      }

      civ_plus = refgrad_plus->civectors()->copy();
      civ_plus->rotate (task_->msrot());
      civ_plus->match(civ_ref);

      displ->element(j,i) = -2.0 * dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();

      refgrad_minus = make_shared<Reference> (*ref_, nullptr);
      refgrad_minus = nullptr;
      
      task_ = std::make_shared<CASPT2Ener>(idata_out, geom_, refgrad_minus);
      task_->compute();
      refgrad_minus  = task_->conv_to_ref();

      acoeff_minus = make_shared<Matrix>(task_->coeff()->slice(nclosed, nocc));
      for (int im = 0; im != acoeff_ref->mdim(); ++im) {
        double dmatch = blas::dot_product(acoeff_ref->element_ptr(0, im), acoeff_ref->ndim(), acoeff_minus->element_ptr(0,im));
        if (dmatch < 0.0) 
          blas::scale_n(-1.0, acoeff_minus->element_ptr(0, im), acoeff_ref->ndim());
      }

      civ_minus = refgrad_minus->civectors()->copy();
      civ_minus->rotate (task_->msrot());
      civ_minus->match(civ_ref);

      civ_diff = civ_plus->copy();
      *civ_diff -= *civ_minus;
      civ_diff->scale(1.0 / (2.0 * dx_));
      civ_diff->print(/*sort=*/false);
      acoeff_diff = make_shared<Matrix>(*acoeff_plus - *acoeff_minus);
      acoeff_diff->scale(1.0 / (2.0 * dx_));
  
      displ->element(j,i) = dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);

      displ->element(j,i) = 0.0;

      grad->element(j,i) = civ_ref->data(target_state1_)->dot_product(civ_diff->data(target_state2_));
      grad_ci->element(j,i) = grad->element(j,i);

      auto Kfactor = make_shared<Matrix>(*acoeff_ref % *Smn * *acoeff_diff);
      for (int ii = 0; ii != norb; ++ii) {
        for (int ij = 0; ij != norb; ++ij) {
          if (ii != ij) {
            for (auto& iter : civ_ref->det()->phia(ii, ij)) {
              size_t iaA = iter.source;
              size_t iaB = iter.target;
              double sign = static_cast<double>(iter.sign);

              for (size_t ib = 0; ib != lenb; ++ib) {                    
                double factor = civ_ref->data(target_state1_)->data(ib+iaB*lenb) * civ_ref->data(target_state2_)->data(ib+iaA*lenb) * sign;
                grad->element(j,i) += factor * Kfactor->element(ij, ii);
                if ((i + j * 3) == 0)
                  gmo->element(ij, ii) += factor;
              }
            }
            for (size_t ia = 0; ia != lena; ++ia) {
              for (auto& iter : civ_ref->det()->phib(ii, ij)) {
                size_t ibA = iter.source;
                size_t ibB = iter.target;
                double sign = static_cast<double>(iter.sign);
                double factor = civ_ref->data(target_state1_)->data(ibB+ia*lenb) * civ_ref->data(target_state2_)->data(ibA+ia*lenb) * sign;
                grad->element(j,i) += factor * Kfactor->element(ij, ii);
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
  auto grad_basis = make_shared<GradFile>(geom_->natom());
  grad_basis = contract_nacme(nullptr, nullptr, nullptr, nullptr, gfin, /*numerical=*/true);

  grad_ci->print(": CI term", 0);
  grad->print(": All symmetric terms", 0);
  grad_basis->print(": Basis set derivative (analytically calculated)", 0);

  *grad += *grad_basis;
  grad->print(": NACME calculated with finite difference (only zero - zero)", 0);

  return grad;
#else
  return nullptr;
#endif
}

// end
