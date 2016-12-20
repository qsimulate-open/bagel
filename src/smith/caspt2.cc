//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: caspt2.cc
// Copyright (C) 2016 Toru Shiozaki
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

#include <bagel_config.h>

#include <src/scf/hf/fock.h>
#include <src/grad/cpcasscf.h>
#include <src/grad/gradeval.h>
#include <src/grad/finite.h>
#include <src/multi/casscf/cassecond.h>
#include <src/multi/casscf/casnoopt.h>
#include <src/multi/casscf/qvec.h>
#include <src/smith/smith.h>
#include <src/smith/caspt2grad.h>
#include <src/prop/multipole.h>
#include <src/prop/hyperfine.h>


using namespace std;
using namespace bagel;

CASPT2Energy::CASPT2Energy(shared_ptr<const PTree> inp, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
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
  target_state1_ = inp->get<int>("_target", -1);
  target_state2_ = inp->get<int>("_target2", -1);

  timer.tick_print("Reference calculation");

  cout << endl << "  === DF-CASPT2 calculation ===" << endl << endl;
#else
  throw logic_error("single point CASPT2 require SMITH-generated code. Please compile BAGEL with --enable-smith");
#endif
}


// compute smith and set rdms and ci deriv to a member
void CASPT2Energy::compute() {
#ifdef COMPILE_SMITH
  shared_ptr<PTree> smithinput = idata_->get_child("smith");
  smithinput->put("_grad", false);
  smithinput->put("_nacm", false);
  smithinput->put("_hyperfine", do_hyperfine_);

  if (target_state1_!=-1) {
    smithinput->put("_target", target_state1_);
    smithinput->put("_target2", target_state2_);
  }
  auto smith = make_shared<Smith>(smithinput, ref_->geom(), ref_);
  smith->compute();
  
  energy_  = smith->algo()->energyvec();
  ref_->energy() = energy_;

  if (target_state1_!=-1) {
    ncore_   = smith->algo()->info()->ncore();
    coeff_   = smith->coeff();
    msrot_   = smith->msrot();
    xmsrot_  = smith->xmsrot();
    heffrot_ = smith->heffrot();

    auto d1set = [this](shared_ptr<const Matrix> d1t) {
      if (!ncore_) {
        return d1t->copy();
      } else {
        auto out = make_shared<Matrix>(coeff_->mdim(), coeff_->mdim());
        out->copy_block(ncore_, ncore_, coeff_->mdim()-ncore_, coeff_->mdim()-ncore_, d1t);
        return out;
      }
    };
    msrot_->print("msrot = ", 20);
    auto vd1tmp = make_shared<Matrix>(*smith->vd1());
    vd1_ = d1set(vd1tmp);
  }
#endif
}

template<>
shared_ptr<GradFile> FiniteNacm<CASPT2Energy>::compute() {
#ifdef COMPILE_SMITH
  cout << "  NACME evaluation with respect to " << geom_->natom() * 3 << " DOFs" << endl;
  cout << "  Finite difference size (dx) is " << setprecision(8) << dx_ << " Bohr(s)" << endl;

#if 1
  auto displ = make_shared<XYZFile>(geom_->natom());
  auto grad_ci = std::make_shared<GradFile>(geom_->natom());
  auto grad_csf = std::make_shared<GradFile>(geom_->natom());
  auto grad = make_shared<GradFile>(geom_->natom());
  auto vd1a = make_shared<Matrix>(*task_->vd1());

  displ->scale(0.0);
  
  shared_ptr<Dvec> civ_ref = ref_->civectors()->copy();
  int nclosed = ref_->nclosed();
  int nocc = ref_->nocc();
  shared_ptr<const Matrix> acoeff_ref, coeff_ref;
  acoeff_ref = make_shared<Matrix>(task_->coeff()->slice(nclosed, nocc));
  coeff_ref  = make_shared<Matrix>(*task_->coeff());

  civ_ref->rotate (task_->msrot());
  civ_ref->print (/*sort=*/false);
  
  shared_ptr<const Reference> refgrad_plus;
  shared_ptr<const Reference> refgrad_minus;
  shared_ptr<Dvec> civ_plus;
  shared_ptr<Dvec> civ_minus;
  shared_ptr<Dvec> civ_diff;
  shared_ptr<Matrix> acoeff_diff;
  shared_ptr<Matrix> coeff_plus;
  shared_ptr<Matrix> coeff_minus;
  shared_ptr<Matrix> coeff_diff;

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
      mute_stdcout();
      displ->element(j,i) = dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();

      refgrad_plus = make_shared<Reference>(*ref_, nullptr);
      refgrad_plus = refgrad_plus->project_coeff(geom_);

      task_ = std::make_shared<CASPT2Energy>(idata_out, geom_, refgrad_plus);
      task_->compute();
      refgrad_plus  = task_->conv_to_ref();

      coeff_plus = make_shared<Matrix>(*task_->coeff());

      for (int im = 0; im != coeff_ref->mdim(); ++im) {
        double dmatch = blas::dot_product(coeff_ref->element_ptr(0, im), coeff_ref->ndim(), coeff_plus->element_ptr(0,im));
        if (dmatch < 0.0) {
          blas::scale_n(-1.0, coeff_plus->element_ptr(0, im), coeff_ref->ndim());
        }
      }

      civ_plus = refgrad_plus->civectors()->copy();
      civ_plus->rotate(task_->msrot());
      civ_plus->match(civ_ref);

      displ->element(j,i) = -2.0 * dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();

      refgrad_minus = make_shared<Reference>(*ref_, nullptr);
      refgrad_minus = refgrad_minus->project_coeff(geom_);
      
      task_ = std::make_shared<CASPT2Energy>(idata_out, geom_, refgrad_minus);
      task_->compute();
      refgrad_minus  = task_->conv_to_ref();

      coeff_minus = make_shared<Matrix>(*task_->coeff());

      for (int im = 0; im != coeff_ref->mdim(); ++im) {
        double dmatch = blas::dot_product(coeff_ref->element_ptr(0, im), coeff_ref->ndim(), coeff_minus->element_ptr(0,im));
        if (dmatch < 0.0) {
          blas::scale_n(-1.0, coeff_minus->element_ptr(0, im), coeff_ref->ndim());
        }
      }

      civ_minus = refgrad_minus->civectors()->copy();
      civ_minus->rotate(task_->msrot());
      civ_minus->match(civ_ref);

      civ_diff = civ_plus->copy();
      *civ_diff -= *civ_minus;
      civ_diff->scale(1.0 / (2.0 * dx_));
      civ_diff->print(/*sort=*/false);
      coeff_diff = make_shared<Matrix>(*coeff_plus - *coeff_minus);
      coeff_diff->scale(1.0 / (2.0 * dx_));
      acoeff_diff = make_shared<Matrix>(coeff_diff->slice(nclosed, nocc));
  
      displ->element(j,i) = dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);

      displ->element(j,i) = 0.0;

      grad->element(j,i) = civ_ref->data(target_state1_)->dot_product(civ_diff->data(target_state2_));
      grad_ci->element(j,i) = grad->element(j,i);
      grad->element(j,i) = 0.0;
      grad_csf->element(j,i) = 0.0;

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
//                grad->element(j,i) += factor * Uij->element(ij, ii);
                if ((i + j * 3) == 0) {
                  gmo->element(ij, ii) += factor * .5;
                  gmo->element(ii, ij) -= factor * .5;
//                  gmo->element(ij, ii) += factor;
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
//                grad->element(j,i) += factor * Uij->element(ij, ii); 
                if ((i + j * 3) == 0) {
                  gmo->element(ij, ii) += factor * .5;
                  gmo->element(ii, ij) -= factor * .5;
//                  gmo->element(ij, ii) += factor;
                }
              }
            }
          }
        }
      }
    
      const int nmobasis = task_->coeff()->ndim();
      auto Ifactor = make_shared<Matrix>(*coeff_ref % *Smn * *coeff_diff);
      for (int ii = 0; ii != nmobasis; ++ii) {
        for (int ij = 0; ij != nmobasis; ++ij) {
          grad_csf->element(j,i) += vd1a->element(ij, ii) * Ifactor->element(ij, ii);
        }
      }
      resume_stdcout();
      cout << "Finite difference evaluation " << setw(5) << i*3+j+1 << " / " << geom_->natom() * 3 << endl;
    }
  }


  auto gfin = make_shared<Matrix>((*acoeff_ref * (*gmo) ^ *acoeff_ref) + (*coeff_ref * (*vd1a) ^ *coeff_ref));
  auto grad_basis = make_shared<GradFile>(geom_->natom());
  grad_basis = contract_nacme(nullptr, nullptr, nullptr, nullptr, gfin, /*numerical=*/true);

//  grad_ci->print(": CI term, <cJ | (cI(X+dX) - cI(X-dX))>", 0);
//  grad->print(": CI term, Orbital-dependent", 0);
  *grad += *grad_ci;
  grad->print(": CI term", 0);
  grad_csf->print(": First-order CSF term", 0);
  *grad += *grad_csf;
//  grad->print(": CI term + First-order CSF term", 0);
  grad_basis->print(": Basis set derivative (analytically calculated)", 0);

  *grad += *grad_basis;
  grad->print(": NACME calculated with finite difference", 0);
  return grad;
#else
  auto displ = make_shared<XYZFile>(geom_->natom());
  auto grad_msrot = std::make_shared<GradFile>(geom_->natom());
  auto grad_fock = std::make_shared<GradFile>(geom_->natom());
  auto grad_ci = std::make_shared<GradFile>(geom_->natom());
  auto grad_csf = std::make_shared<GradFile>(geom_->natom());
  auto grad = make_shared<GradFile>(geom_->natom());
  auto vd1a = make_shared<Matrix>(*task_->vd1());

  displ->scale(0.0);
  
  shared_ptr<Dvec> civ_ref = ref_->civectors()->copy();
  int nclosed = ref_->nclosed();
  int nocc = ref_->nocc();
  shared_ptr<const Matrix> acoeff_ref, coeff_ref;
  acoeff_ref = make_shared<Matrix>(task_->coeff()->slice(nclosed, nocc));
  coeff_ref  = make_shared<Matrix>(*task_->coeff());

  auto heffrot_ref = make_shared<Matrix>(*task_->heffrot());
  shared_ptr<Matrix> xmsrot_ref;
  if(task_->xmsrot())
    xmsrot_ref = make_shared<Matrix>(*task_->xmsrot());
  auto msrot_ref   = make_shared<Matrix>(*task_->msrot());

  civ_ref->print (/*sort=*/false);
  
  shared_ptr<const Reference> refgrad_plus;
  shared_ptr<const Reference> refgrad_minus;
  shared_ptr<Dvec> civ_plus;
  shared_ptr<Dvec> civ_minus;
  shared_ptr<Dvec> civ_diff;
  shared_ptr<Matrix> acoeff_diff;
  shared_ptr<Matrix> coeff_plus;
  shared_ptr<Matrix> coeff_minus;
  shared_ptr<Matrix> coeff_diff;
  shared_ptr<Matrix> heffrot_plus;
  shared_ptr<Matrix> heffrot_minus;
  shared_ptr<Matrix> heffrot_diff;
  shared_ptr<Matrix> xmsrot_plus;
  shared_ptr<Matrix> xmsrot_minus;
  shared_ptr<Matrix> xmsrot_diff;
//  std::vector<std::shared_ptr<MultiTensor>> t2all_plus;
//  std::vector<std::shared_ptr<MultiTensor>> t2all_minus;
//  std::vector<std::shared_ptr<MultiTensor>> t2all_diff;

  auto Smn = make_shared<Overlap>(geom_);
  const int norb = civ_ref->det()->norb();
  const int nmobasis = task_->coeff()->ndim();
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

      refgrad_plus = make_shared<Reference>(*ref_, nullptr);
      refgrad_plus = refgrad_plus->project_coeff(geom_);

      task_ = std::make_shared<CASPT2Energy>(idata_out, geom_, refgrad_plus);
      task_->compute();
      refgrad_plus  = task_->conv_to_ref();

      coeff_plus = make_shared<Matrix>(*task_->coeff());
      heffrot_plus = make_shared<Matrix>(*task_->heffrot());
      if(task_->xmsrot())
        xmsrot_plus  = make_shared<Matrix>(*task_->xmsrot());

      for (int im = 0; im != coeff_ref->mdim(); ++im) {
        double dmatch = blas::dot_product(coeff_ref->element_ptr(0, im), coeff_ref->ndim(), coeff_plus->element_ptr(0,im));
        if (dmatch < 0.0) {
          blas::scale_n(-1.0, coeff_plus->element_ptr(0, im), coeff_ref->ndim());
        }
      }

      for (int im = 0; im != heffrot_ref->mdim(); ++im) {
        for (int in = 0; in != heffrot_ref->ndim(); ++in) {
          double dmatch = heffrot_ref->element(im, in) * heffrot_plus->element(im, in);
          if (dmatch < 0.0) heffrot_plus->element(im, in) *= -1.0;
        }
      }
      
      if(task_->xmsrot())
        for (int im = 0; im != xmsrot_ref->mdim(); ++im) {
          for (int in = 0; in != xmsrot_ref->ndim(); ++in) {
            double dmatch = xmsrot_ref->element(im, in) * xmsrot_plus->element(im, in);
            if (dmatch < 0.0) xmsrot_plus->element(im, in) *= -1.0;
          }
        }


      civ_plus = refgrad_plus->civectors()->copy();
      civ_plus->match(civ_ref);

      displ->element(j,i) = -2.0 * dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);
      geom_->print_atoms();

      refgrad_minus = make_shared<Reference>(*ref_, nullptr);
      refgrad_minus = refgrad_minus->project_coeff(geom_);
      
      task_ = std::make_shared<CASPT2Energy>(idata_out, geom_, refgrad_minus);
      task_->compute();
      refgrad_minus  = task_->conv_to_ref();

      coeff_minus = make_shared<Matrix>(*task_->coeff());
      heffrot_minus = make_shared<Matrix>(*task_->heffrot());
      if(task_->xmsrot())
        xmsrot_minus  = make_shared<Matrix>(*task_->xmsrot());

      for (int im = 0; im != coeff_ref->mdim(); ++im) {
        double dmatch = blas::dot_product(coeff_ref->element_ptr(0, im), coeff_ref->ndim(), coeff_minus->element_ptr(0,im));
        if (dmatch < 0.0) {
          blas::scale_n(-1.0, coeff_minus->element_ptr(0, im), coeff_ref->ndim());
        }
      }

      for (int im = 0; im != heffrot_ref->mdim(); ++im) {
        for (int in = 0; in != heffrot_ref->ndim(); ++in) {
          double dmatch = heffrot_ref->element(im, in) * heffrot_minus->element(im, in);
          if (dmatch < 0.0) heffrot_minus->element(im, in) *= -1.0;
        }
      }

      if(task_->xmsrot())
        for (int im = 0; im != xmsrot_ref->mdim(); ++im) {
          for (int in = 0; in != xmsrot_ref->ndim(); ++in) {
            double dmatch = xmsrot_ref->element(im, in) * xmsrot_minus->element(im, in);
            if (dmatch < 0.0) xmsrot_minus->element(im, in) *= -1.0;
          }
        }

      civ_minus = refgrad_minus->civectors()->copy();
      civ_minus->match(civ_ref);

      civ_diff = civ_plus->copy();
      *civ_diff -= *civ_minus;
      civ_diff->scale(1.0 / (2.0 * dx_));
      civ_diff->print(/*sort=*/false);
      coeff_diff = make_shared<Matrix>(*coeff_plus - *coeff_minus);
      coeff_diff->scale(1.0 / (2.0 * dx_));
      acoeff_diff = make_shared<Matrix>(coeff_diff->slice(nclosed, nocc));
      heffrot_diff = make_shared<Matrix>(*heffrot_plus - *heffrot_minus);
      if(task_->xmsrot())
        xmsrot_diff  = make_shared<Matrix>(*xmsrot_plus - *xmsrot_minus);
      heffrot_diff->scale(1.0 / (2.0 * dx_));
      if(task_->xmsrot())
        xmsrot_diff->scale(1.0 / (2.0 * dx_));
  
      displ->element(j,i) = dx_;
      geom_ = make_shared<Geometry>(*geom_, displ);

      displ->element(j,i) = 0.0;

      grad_msrot->element(j,i) = blas::dot_product(heffrot_ref->element_ptr(0, target_state1_), heffrot_ref->ndim(), heffrot_diff->element_ptr(0, target_state2_));
      grad_fock->element(j,i) = 0.0;
      grad_ci->element(j,i) = 0.0;
      grad_csf->element(j,i) = 0.0;

      // Fock term
      if(task_->xmsrot())
        for (int im = 0; im != heffrot_ref->mdim(); ++im)
          for (int in = 0; in != heffrot_ref->ndim(); ++in) {
            if (im==in) continue;
            else {
              grad_fock->element(j,i) += blas::dot_product(xmsrot_ref->element_ptr(0, im),
                  heffrot_ref->ndim(), xmsrot_diff->element_ptr(0, in)) *
                (heffrot_ref->element(im, target_state1_) * heffrot_ref->element(in, target_state2_));
            }
          }

      // CI term
      for (int im = 0; im != msrot_ref->mdim(); ++im) {
        for (int in = 0; in != msrot_ref->ndim(); ++in) {
          if (im==in) continue;
          else {
            grad_ci->element(j,i) += civ_ref->data(im)->dot_product(civ_diff->data(in)) *
              (msrot_ref->element(im, target_state1_) * msrot_ref->element(in, target_state2_));
          }
        }
      }

      auto Uij = make_shared<Matrix>(*acoeff_ref % *Smn * *acoeff_diff);
      Uij->print("Uij = ", 20);
      for (int ii = 0; ii != norb; ++ii) {
        for (int ij = 0; ij != norb; ++ij) {
          if (ii != ij) {
            for (auto& iter : civ_ref->det()->phia(ii, ij)) {
              size_t iaA = iter.source;
              size_t iaB = iter.target;
              double sign = static_cast<double>(iter.sign);

              for (size_t ib = 0; ib != lenb; ++ib) {
                for (int im = 0; im != msrot_ref->mdim(); ++im) {
                  for (int in = 0; in != msrot_ref->ndim(); ++in) {
                    if (im==in) continue;
                    else {
                      double factor = civ_ref->data(im)->data(ib+iaB*lenb) * civ_ref->data(in)->data(ib+iaA*lenb) * sign * 
                                    (msrot_ref->element(im, target_state1_) * msrot_ref->element(in, target_state2_)
                                    -msrot_ref->element(in, target_state1_) * msrot_ref->element(im, target_state2_)) * .5;
                      grad_ci->element(j,i) += factor * (Uij->element(ij, ii) - Uij->element(ii, ij)) * .5;
                      if ((i + j * 3) == 0) {
                        gmo->element(ij, ii) += factor * .5;
                        gmo->element(ii, ij) -= factor * .5;
                      }
                    }
                  }
                }
              }
            }
            for (size_t ia = 0; ia != lena; ++ia) {
              for (auto& iter : civ_ref->det()->phib(ii, ij)) {
                size_t ibA = iter.source;
                size_t ibB = iter.target;
                double sign = static_cast<double>(iter.sign);
                for (int im = 0; im != msrot_ref->mdim(); ++im) {
                  for (int in = 0; in != msrot_ref->ndim(); ++in) {
                    if (im==in) continue;
                    else {
                      double factor = civ_ref->data(im)->data(ibB+ia*lenb) * civ_ref->data(in)->data(ibA+ia*lenb) * sign * 
                                    (msrot_ref->element(im, target_state1_) * msrot_ref->element(in, target_state2_)
                                    -msrot_ref->element(in, target_state1_) * msrot_ref->element(im, target_state2_)) * .5;
                      grad_ci->element(j,i) += factor * (Uij->element(ij, ii) - Uij->element(ii, ij)) * .5;
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
      }
      auto Ifactor = make_shared<Matrix>(*coeff_ref % *Smn * *coeff_diff);
      // additional part due to perturbation
      for (int ii = 0; ii != nmobasis; ++ii) {
        for (int ij = 0; ij != nmobasis; ++ij) {
          if (ii != ij) {
            grad_csf->element(j,i) += vd1a->element(ij, ii) * Ifactor->element(ij, ii);
          }
        }
      }
    }
  }

  auto gfin = make_shared<Matrix>((*acoeff_ref * (*gmo) ^ *acoeff_ref) + (*coeff_ref * (*vd1a) ^ *coeff_ref));
  auto grad_basis = make_shared<GradFile>(geom_->natom());
  grad_basis = contract_nacme(nullptr, nullptr, nullptr, nullptr, gfin, /*numerical=*/true);

  grad_msrot->print(": MSRot term", 0);
  grad_ci->print(": CI term", 0);
  grad_fock->print(": Fock term", 0);
  *grad_msrot += *grad_ci;
  *grad_msrot += *grad_fock;
  *grad += *grad_msrot;
  grad->print(": MSROT + CI + Fock", 0);
  grad_csf->print(": First Order CSF term", 0);
  *grad += *grad_csf;
  grad_basis->print(": Basis set derivative (analytically calculated)", 0);

  *grad += *grad_basis;
  grad->print(": NACME calculated with finite difference", 0);

  return grad;
#endif
#else
  return nullptr;
#endif
}

shared_ptr<const Reference> CASPT2Energy::conv_to_ref() const {
 return std::make_shared<Reference>(ref_->geom(), ref_->coeff(), ref_->nclosed(), ref_->nact(), ref_->nvirt(), energy_,
                               fci_->rdm1(), fci_->rdm2(), fci_->rdm1_av(), fci_->rdm2_av(), fci_->conv_to_ciwfn());
}
