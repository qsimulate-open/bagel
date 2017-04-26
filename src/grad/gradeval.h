//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gradeval.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_GRAD_GRADEVAL_H
#define __SRC_GRAD_GRADEVAL_H

#include <src/scf/hf/rohf.h>
#include <src/scf/ks/ks.h>
#include <src/scf/dhf/dirac.h>
#include <src/multi/casscf/casscf.h>
#include <src/pt2/mp2/mp2grad.h>
#include <src/grad/gradeval_base.h>
#include <src/smith/caspt2grad.h>

// T should have
// o Constructor with the input and geometry
// o void compute()
// o std::shared_ptr<Reference> conv_to_ref()

namespace bagel {

template<typename T>
class GradEval : public GradEval_base {
  protected:
    std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Reference> ref_;

    std::shared_ptr<T> task_;

    double energy_;
    int target_state_;

    int maxziter_;

    void init() {
      if (geom_->external())
        throw std::logic_error("Gradients with external fields have not been implemented.");
      // target has to be passed to T (for CASPT2, for instance)
      auto idata_out = std::make_shared<PTree>(*idata_);
      idata_out->put("_target", target_state_);
      idata_out->put("_maxziter", maxziter_);
      task_ = std::make_shared<T>(idata_out, geom_, ref_);
      task_->compute();
      ref_  = task_->conv_to_ref();
      energy_ = ref_->energy(target_state_);
      geom_ = ref_->geom();
    }

  public:
    // Constructor performs energy calculation
    GradEval(std::shared_ptr<const PTree> idata, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref, const int target, const int maxziter = 100)
      : GradEval_base(geom), idata_(idata), ref_(ref), target_state_(target), maxziter_(maxziter) {
      init();
    }

    // compute() computes effective density matrices and perform gradient contractions
    std::shared_ptr<GradFile> compute() { throw std::logic_error("Nuclear gradient for this method has not been implemented"); }

    double energy() const { return energy_; }

    std::shared_ptr<const Reference> ref() const { return ref_; }
};

// specialization
template<> std::shared_ptr<GradFile> GradEval<RHF>::compute();
template<> std::shared_ptr<GradFile> GradEval<UHF>::compute();
template<> std::shared_ptr<GradFile> GradEval<ROHF>::compute();
template<> std::shared_ptr<GradFile> GradEval<KS>::compute();
template<> std::shared_ptr<GradFile> GradEval<MP2Grad>::compute();
template<> std::shared_ptr<GradFile> GradEval<Dirac>::compute();
template<> std::shared_ptr<GradFile> GradEval<CASPT2Grad>::compute();

// CASSCF is slightly more complicated. These functions are implemented in casgrad.cc
template<> void GradEval<CASSCF>::init();
template<> std::shared_ptr<GradFile> GradEval<CASSCF>::compute();

template<typename T>
class NacmEval : public GradEval_base {
  protected:
    std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Reference> ref_;

    std::shared_ptr<T> task_;

    double energy1_;
    double energy2_;
    int target_state1_;
    int target_state2_;
    // nacmtype is 0 (full derivative coupling), 1 (interstate coupling) or 2 (full derivative coupling with built-in ETF factor: see Subotnik)
    int nacmtype_;

    int maxziter_;

    void init() {
      if (geom_->external())
        throw std::logic_error("Gradients with external fields have not been implemented.");
      auto idata_out = std::make_shared<PTree>(*idata_);
      idata_out->put("_target", target_state1_);
      idata_out->put("_target2", target_state2_);
      idata_out->put("_nacmtype", nacmtype_);
      idata_out->put("_maxziter", maxziter_);
      task_ = std::make_shared<T>(idata_out, geom_, ref_);
      task_->compute();
      ref_  = task_->conv_to_ref();
      energy1_ = ref_->energy(target_state1_);
      energy2_ = ref_->energy(target_state2_);
      geom_ = ref_->geom();
    }

  public:
    NacmEval(std::shared_ptr<const PTree> idata, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref, const int target1, const int target2, const int nacmtype, const int maxziter = 100)
      : GradEval_base(geom), idata_(idata), ref_(ref), target_state1_(target1), target_state2_(target2), nacmtype_(nacmtype), maxziter_(maxziter) {
      init();
    }

    // compute() computes effective density matrices and perform gradient contractions
    std::shared_ptr<GradFile> compute() { throw std::logic_error("NACME for this method has not been implemented"); }

    double energy1 () const { return energy1_; }
    double energy2 () const { return energy2_; }

    std::shared_ptr<const Reference> ref() const { return ref_; }
};
template<> std::shared_ptr<GradFile> NacmEval<CASPT2Nacm>::compute();

// CASSCF code for NACME is basically same to the gradient one, but little different due to some additional terms...
template<> void NacmEval<CASSCF>::init();
template<> std::shared_ptr<GradFile> NacmEval<CASSCF>::compute();

template<typename T>
class DgradEval : public GradEval_base {
  protected:
    std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Reference> ref_;

    std::shared_ptr<T> task_;

    double energy1_;
    double energy2_;
    int target_state1_;
    int target_state2_;

    int maxziter_;

    void init() {
      if (geom_->external())
        throw std::logic_error("Gradients with external fields have not been implemented.");
      auto idata_out = std::make_shared<PTree>(*idata_);
      idata_out->put("_target", target_state1_);
      idata_out->put("_target2", target_state2_);
      idata_out->put("_maxziter", maxziter_);
      task_ = std::make_shared<T>(idata_out, geom_, ref_);
      task_->compute();
      ref_  = task_->conv_to_ref();
      energy1_ = ref_->energy(target_state1_);
      energy2_ = ref_->energy(target_state2_);
      geom_ = ref_->geom();
    }

  public:
    DgradEval(std::shared_ptr<const PTree> idata, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref, const int target1, const int target2, const int maxziter)
      : GradEval_base(geom), idata_(idata), ref_(ref), target_state1_(target1), target_state2_(target2), maxziter_(maxziter) {
      init();
    }

    // compute() computes effective density matrices and perform gradient contractions
    std::shared_ptr<GradFile> compute() { throw std::logic_error("State difference gradient for this method has not been implemented"); }

    double energy1 () const { return energy1_; }
    double energy2 () const { return energy2_; }

    std::shared_ptr<const Reference> ref() const { return ref_; }
};
template<> void DgradEval<CASSCF>::init();
template<> std::shared_ptr<GradFile> DgradEval<CASSCF>::compute();


}

#endif
