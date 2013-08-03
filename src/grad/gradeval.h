//
// BAGEL - Parallel electron correlation program.
// Filename: gradeval.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_GRAD_GRADEVAL_H
#define __SRC_GRAD_GRADEVAL_H

#include <src/scf/rohf.h>
#include <src/ks/ks.h>
#include <src/mp2/mp2grad.h>
#include <src/grad/gradeval_base.h>
#include <src/casscf/werner.h>
#include <src/casscf/supercigrad.h>
#include <src/rel/dirac.h>
#include <src/rel/dmp2grad.h>

// T should have
// o Constructor with the input and geometry
// o void compute()
// o std::shared_ptr<Reference> conv_to_ref()

namespace bagel {

template<typename T>
class GradEval : public GradEval_base {
  protected:
    std::shared_ptr<const Reference> ref_;

    std::shared_ptr<T> task_;

    double energy_;

  public:
    // Constructor performs energy calculation
    GradEval(std::shared_ptr<const PTree> idata, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref) : GradEval_base(geom) {
      if (geom->external()) throw std::logic_error("Gradients with external fields have not been implemented.");
      task_ = std::make_shared<T>(idata, geom, ref);
      task_->compute();
      ref_  = task_->conv_to_ref();
      energy_ = ref_->energy();
      geom_ = ref_->geom();
    }

    // compute() computes effective density matrices and perform gradient contractions
    std::shared_ptr<GradFile> compute() { throw std::logic_error("Nuclear gradient for this method has not been implemented"); }

    double energy() const { return energy_; }

    std::shared_ptr<const Reference> ref() const { return ref_; }
};

// specialization
template<> std::shared_ptr<GradFile> GradEval<SCF>::compute();
template<> std::shared_ptr<GradFile> GradEval<UHF>::compute();
template<> std::shared_ptr<GradFile> GradEval<ROHF>::compute();
template<> std::shared_ptr<GradFile> GradEval<KS>::compute();
template<> std::shared_ptr<GradFile> GradEval<MP2Grad>::compute();
template<> std::shared_ptr<GradFile> GradEval<WernerKnowles>::compute();
template<> std::shared_ptr<GradFile> GradEval<SuperCI>::compute();
template<> std::shared_ptr<GradFile> GradEval<SuperCIGrad>::compute();
template<> std::shared_ptr<GradFile> GradEval<Dirac>::compute();
template<> std::shared_ptr<GradFile> GradEval<DMP2Grad>::compute();

}

#endif
