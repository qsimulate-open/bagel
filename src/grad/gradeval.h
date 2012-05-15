//
// Newint - Parallel electron correlation program.
// Filename: gradeval.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_GRAD_GRADEVAL_H
#define __SRC_GRAD_GRADEVAL_H

#include <map>
#include <string>
#include <vector>
#include <memory>
#include <src/wfn/reference.h>
#include <src/scf/scf.h>
#include <src/scf/uhf.h>
#include <src/grad/gradeval_base.h>

// T should have
// o Constructor with the input and geometry
// o void compute()
// o std::shared_ptr<Referenc> conv_to_ref()

template<typename T>
class GradEval : public GradEval_base {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<Reference> ref_;

    double energy_;

  public:
    // Constructor performs energy calculation
    GradEval(std::multimap<std::string, std::string>& idata, const std::shared_ptr<const Geometry> geom) : GradEval_base(geom), geom_(geom) {
      T task(idata, geom);
      task.compute();
      ref_  = task.conv_to_ref();
      energy_ = ref_->energy();
    };
    ~GradEval() {};

    // compute() computes effective density matrices and perform gradient contractions
    std::shared_ptr<GradFile> compute() const { assert(false); };

    double energy() const { return energy_; };
}; 

// specialization
template<> std::shared_ptr<GradFile> GradEval<SCF<1> >::compute() const;
template<> std::shared_ptr<GradFile> GradEval<UHF>::compute() const;

#endif
