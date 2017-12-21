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
#include <src/grad/gradinfo.h>
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
    std::vector<double> dipole_;

    void init() {
      if (geom_->external())
        throw std::logic_error("Gradients with external fields have not been implemented.");
      auto idata_out = std::make_shared<PTree>(*idata_);
      task_ = std::make_shared<T>(idata_out, geom_, ref_);
      task_->compute();
      ref_  = task_->conv_to_ref();
      geom_ = ref_->geom();
    }

  public:

    // Constructor performs energy calculation
    GradEval(std::shared_ptr<const PTree> idata, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref)
      : GradEval_base(geom), idata_(idata), ref_(ref) {
      init();
    }

    // computes effective density matrices and perform gradient contractions
    std::shared_ptr<GradFile> compute(const std::string jobtitle = "force", std::shared_ptr<const GradInfo> gradinfo = std::make_shared<const GradInfo>())
      { throw std::logic_error("Nuclear gradient for this method has not been implemented"); }

    // computes unrelaxed dipole moments
    void compute_dipole() const
      { throw std::logic_error("compute_dipole() only works for CASSCF and CASPT2"); }

    double energy() const { return energy_; }
    std::vector<double> energyvec() const { return ref_->energy(); }
    const std::vector<double>& dipole() const { return dipole_; }
    double dipole(int i) const { return dipole_[i]; }

    std::shared_ptr<const Reference> ref() const { return ref_; }
};

// specialization
template<> std::vector<double> GradEval<MP2Grad>::energyvec() const;
template<> std::vector<double> GradEval<CASPT2Grad>::energyvec() const;

template<> std::shared_ptr<GradFile> GradEval<RHF>::compute(const std::string jobtitle, std::shared_ptr<const GradInfo> gradinfo);
template<> std::shared_ptr<GradFile> GradEval<UHF>::compute(const std::string jobtitle, std::shared_ptr<const GradInfo> gradinfo);
template<> std::shared_ptr<GradFile> GradEval<ROHF>::compute(const std::string jobtitle, std::shared_ptr<const GradInfo> gradinfo);
template<> std::shared_ptr<GradFile> GradEval<KS>::compute(const std::string jobtitle, std::shared_ptr<const GradInfo> gradinfo);
template<> std::shared_ptr<GradFile> GradEval<MP2Grad>::compute(const std::string jobtitle, std::shared_ptr<const GradInfo> gradinfo);
template<> std::shared_ptr<GradFile> GradEval<Dirac>::compute(const std::string jobtitle, std::shared_ptr<const GradInfo> gradinfo);
template<> std::shared_ptr<GradFile> GradEval<CASPT2Grad>::compute(const std::string jobtitle, std::shared_ptr<const GradInfo> gradinfo);

// These functions are implemented in casgrad.cc
template<> void GradEval<CASSCF>::init();
template<> std::shared_ptr<GradFile> GradEval<CASSCF>::compute(const std::string jobtitle, std::shared_ptr<const GradInfo> gradinfo);

template<> void GradEval<CASSCF>::compute_dipole() const;
template<> void GradEval<CASPT2Grad>::compute_dipole() const;

}

#endif
