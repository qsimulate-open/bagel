//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: finite.h
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


#ifndef __SRC_GRAD_FINITE_H
#define __SRC_GRAD_FINITE_H

#include <src/multi/casscf/casscf.h>
#include <src/smith/caspt2grad.h>
#include <src/grad/gradeval_base.h>

namespace bagel {

class FiniteGrad : public GradEval_base {
  protected:
    std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Reference> ref_;

    mutable std::shared_ptr<Muffle> muffle_;
    double energy_;

    int target_state_;
    double dx_;
    int nproc_;

  public:
    // Constructor does nothing here
    FiniteGrad(std::shared_ptr<const PTree> idata, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref, const int target, const double dx, const int nproc)
      : GradEval_base(geom), idata_(idata), ref_(ref), target_state_(target), dx_(dx), nproc_(nproc) {
    }

    std::shared_ptr<GradFile> compute();

    double energy() const { return energy_; }

    std::shared_ptr<const Reference> ref() const { return ref_; }
};

template<typename T>
class FiniteNacm : public GradEval_base {
  protected:
    std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Reference> ref_;

    mutable std::shared_ptr<Muffle> muffle_;

    std::shared_ptr<T> task_;

    double energy1_;
    double energy2_;

    int target_state1_;
    int target_state2_;
    double dx_;
    int nproc_;

    void init() {
      if (geom_->external())
        throw std::logic_error("Gradients with external fields have not been implemented.");
      auto idata_out = std::make_shared<PTree>(*idata_);
      idata_out->put("_target", target_state1_);
      idata_out->put("_target2", target_state2_);
      task_ = std::make_shared<T>(idata_out, geom_, ref_);
      task_->compute();
      ref_  = task_->conv_to_ref();
      energy1_ = task_->energy(target_state1_);
      energy2_ = task_->energy(target_state2_);
      std::cout << std::setprecision(8) << "  Energy = " << energy1_ << " and " << energy2_ << std::endl;
      geom_ = ref_->geom();
    }

  public:
    FiniteNacm(std::shared_ptr<const PTree> idata, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref, const int target, const int target2, const double dx, const int nproc)
      : GradEval_base(geom), idata_(idata), ref_(ref), target_state1_(target), target_state2_(target2), dx_(dx), nproc_(nproc) {
      init();
    }

    std::shared_ptr<GradFile> compute() { throw std::logic_error("NACME for this method has not been implemented"); }

    double energy1 () const { return energy1_; }
    double energy2 () const { return energy2_; }

    std::shared_ptr<const Reference> ref() const { return ref_; }
};

template<> void FiniteNacm<CASSCF>::init();
template<> std::shared_ptr<GradFile> FiniteNacm<CASSCF>::compute();
template<> std::shared_ptr<GradFile> FiniteNacm<CASPT2Energy>::compute();

}
#endif
