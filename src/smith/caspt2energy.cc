//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: caspt2energy.cc
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
    smith->compute_gradient(target_state1_, target_state2_, make_shared<const NacmType>());
    ncore_   = smith->algo()->info()->ncore();
    coeff_   = smith->coeff();
    msrot_   = smith->msrot();

    auto d1set = [this](shared_ptr<const Matrix> d1t) {
      if (!ncore_) {
        return d1t->copy();
      } else {
        auto out = make_shared<Matrix>(coeff_->mdim(), coeff_->mdim());
        out->copy_block(ncore_, ncore_, coeff_->mdim()-ncore_, coeff_->mdim()-ncore_, d1t);
        return out;
      }
    };
    auto vd1tmp = make_shared<Matrix>(*smith->vd1());
    vd1_ = d1set(vd1tmp);
  }
#endif
}

template<>
shared_ptr<GradFile> FiniteNacm<CASPT2Energy>::compute() {
#ifdef COMPILE_SMITH
  const int natom = geom_->natom();
  cout << "  Derivative coupling evaluation with respect to " << natom * 3 << " DOFs" << endl;
  cout << "  Finite difference size (dx) is " << setprecision(8) << dx_ << " Bohr" << endl;

  shared_ptr<Dvec> civ_ref = ref_->civectors()->copy();
  civ_ref->rotate (task_->msrot());
  civ_ref->print (/*sort=*/false);

  const int nclosed = ref_->nclosed();
  const int nocc = ref_->nocc();

  Timer timer;
  muffle_ = make_shared<Muffle>("finite.log");

  auto grad = make_shared<GradFile>(natom);

  auto acoeff_ref = make_shared<Matrix>(task_->coeff()->slice(nclosed, nocc));
  auto coeff_ref  = make_shared<Matrix>(*task_->coeff());

  const int norb = civ_ref->det()->norb();
  auto gmo = make_shared<Matrix>(norb, norb);
  auto vd1a = make_shared<Matrix>(*task_->vd1());

  auto idata_out = std::make_shared<PTree>(*idata_);
  idata_out->put("_target", target_state1_);
  idata_out->put("_target2", target_state2_);

  const int ncomm = mpi__->world_size() / nproc_;
  const int npass = (natom * 3 - 1) / ncomm + 1;

  for (int ipass = 0; ipass != npass; ++ipass) {
    const int ncolor = (ipass == (npass-1)) ? (natom * 3) % ncomm : ncomm;
    const int icomm = mpi__->world_rank() % ncolor;
    if (ncolor != 0) {
      mpi__->split(ncolor);
    } else {
      continue;
    }

    const int counter = icomm + ncomm * ipass;
    const int i = counter / 3;
    const int j = counter % 3;

    muffle_->mute();
    shared_ptr<Dvec> civ_plus;
    shared_ptr<Matrix> coeff_plus;
    {
      auto displ = make_shared<XYZFile>(natom);
      displ->element(j,i) = dx_;
      auto geom_plus = make_shared<Geometry>(*geom_, displ, idata_, false, false);
      geom_plus->print_atoms();

      shared_ptr<const Reference> ref_plus;
      if (ref_)
        ref_plus = ref_->project_coeff(geom_plus);

      task_ = std::make_shared<CASPT2Energy>(idata_out, geom_plus, ref_plus);
      task_->compute();
      ref_plus  = task_->conv_to_ref();

      coeff_plus = make_shared<Matrix>(*task_->coeff());

      for (int im = 0; im != coeff_ref->mdim(); ++im) {
        double dmatch = blas::dot_product(coeff_ref->element_ptr(0, im), coeff_ref->ndim(), coeff_plus->element_ptr(0,im));
        if (dmatch < 0.0) {
          blas::scale_n(-1.0, coeff_plus->element_ptr(0, im), coeff_ref->ndim());
        }
      }

      civ_plus = ref_plus->civectors()->copy();
      civ_plus->rotate(task_->msrot());
      civ_plus->match(civ_ref);
    }

    shared_ptr<Dvec> civ_minus;
    shared_ptr<Matrix> coeff_minus;
    {
      auto displ = make_shared<XYZFile>(natom);
      displ->element(j,i) = -dx_;
      auto geom_minus = make_shared<Geometry>(*geom_, displ, idata_, false, false);
      geom_minus->print_atoms();

      shared_ptr<const Reference> ref_minus;
      if (ref_)
        ref_minus = ref_->project_coeff(geom_minus);

      task_ = std::make_shared<CASPT2Energy>(idata_out, geom_minus, ref_minus);
      task_->compute();
      ref_minus  = task_->conv_to_ref();

      coeff_minus = make_shared<Matrix>(*task_->coeff());

      for (int im = 0; im != coeff_ref->mdim(); ++im) {
        double dmatch = blas::dot_product(coeff_ref->element_ptr(0, im), coeff_ref->ndim(), coeff_minus->element_ptr(0,im));
        if (dmatch < 0.0) {
          blas::scale_n(-1.0, coeff_minus->element_ptr(0, im), coeff_ref->ndim());
        }
      }

      civ_minus = ref_minus->civectors()->copy();
      civ_minus->rotate(task_->msrot());
      civ_minus->match(civ_ref);
    }

    auto civ_diff = civ_plus->copy();
    *civ_diff -= *civ_minus;
    civ_diff->scale(1.0 / (2.0 * dx_));
    civ_diff->print(/*sort=*/false);
    auto coeff_diff = make_shared<Matrix>(*coeff_plus - *coeff_minus);
    coeff_diff->scale(1.0 / (2.0 * dx_));
    auto acoeff_diff = make_shared<Matrix>(coeff_diff->slice(nclosed, nocc));

    auto Smn = make_shared<Overlap>(geom_);
    auto Uij = make_shared<Matrix>(*acoeff_ref % *Smn * *acoeff_diff);
    if (mpi__->rank() == 0) {
      grad->element(j,i) = civ_ref->data(target_state1_)->dot_product(civ_diff->data(target_state2_));
      for (int ii = 0; ii != norb; ++ii) {
        for (int ij = 0; ij != norb; ++ij) {
          const int lena = civ_ref->det()->lena();
          const int lenb = civ_ref->det()->lenb();
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
      const int nmobasis = task_->coeff()->ndim();
      auto Ifactor = make_shared<Matrix>(*coeff_ref % *Smn * *coeff_diff);
      for (int ii = 0; ii != nmobasis; ++ii) {
        for (int ij = 0; ij != nmobasis; ++ij) {
          grad->element(j,i) += vd1a->element(ij, ii) * Ifactor->element(ij, ii);
        }
      }
    }
    muffle_->unmute();
    mpi__->merge();
    stringstream ss; ss << "Finite difference evaluation (" << setw(2) << i*3+j+1 << " / " << geom_->natom() * 3 << ")";
    timer.tick_print(ss.str());
  }

  grad->allreduce();
  gmo->allreduce();

  auto gfin = make_shared<Matrix>((*acoeff_ref * (*gmo) ^ *acoeff_ref) + (*coeff_ref * (*vd1a) ^ *coeff_ref));
  auto grad_basis = make_shared<GradFile>(geom_->natom());
  grad_basis = contract_gradient(nullptr, nullptr, nullptr, nullptr, gfin, /*numerical=*/true);

  *grad += *grad_basis;
  grad->print(": NACME calculated with finite difference", 0);
  return grad;
#else
  return nullptr;
#endif
}

shared_ptr<const Reference> CASPT2Energy::conv_to_ref() const {
 return std::make_shared<Reference>(ref_->geom(), ref_->coeff(), ref_->nclosed(), ref_->nact(), ref_->nvirt(), energy_,
                               fci_->rdm1(), fci_->rdm2(), fci_->rdm1_av(), fci_->rdm2_av(), fci_->conv_to_ciwfn());
}
