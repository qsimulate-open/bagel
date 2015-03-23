//
// BAGEL - Parallel electron correlation program.
// Filename: spinfreebase.h
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


#ifndef __SRC_SMITH_SPINFREEBASE_H
#define __SRC_SMITH_SPINFREEBASE_H

#include <chrono>
#include <src/ci/fci/civec.h>
#include <src/util/vec.h>
#include <src/smith/denom.h>
#include <src/smith/multitensor.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {

class SpinFreeMethod {
  protected:
    IndexRange virt_;
    IndexRange active_;
    IndexRange closed_;
    IndexRange all_;
    IndexRange ci_;

    // TODO these are redundant.
    std::shared_ptr<const IndexRange> rvirt_;
    std::shared_ptr<const IndexRange> ractive_;
    std::shared_ptr<const IndexRange> rclosed_;
    std::shared_ptr<const IndexRange> rci_;

    std::shared_ptr<const SMITH_Info> ref_;

    std::shared_ptr<const Coeff> coeff_;
    std::shared_ptr<const Civec> civec_;
    double e0_;
    double core_energy_;
    double energy_;

    std::shared_ptr<Tensor> v2_;
    std::shared_ptr<Tensor> f1_;
    std::shared_ptr<Tensor> h1_;

    // contains the current RDMs to be used in smith
    std::shared_ptr<Tensor> rdm1_;
    std::shared_ptr<Tensor> rdm2_;
    std::shared_ptr<Tensor> rdm3_;
    std::shared_ptr<Tensor> rdm4_;

    // contains all the RDMs (for multistate runs)
    std::shared_ptr<Vec<Tensor>> rdm1all_;
    std::shared_ptr<Vec<Tensor>> rdm2all_;
    std::shared_ptr<Vec<Tensor>> rdm3all_;
    std::shared_ptr<Vec<Tensor>> rdm4all_;
    // the function to set RDMs to rdm1_, rdm2_, etc
    void set_rdm(const int jst, const int ist);

    // original determinants (for use in output)
    std::shared_ptr<const Determinants> det_;

    // rdm ci derivatives
    std::shared_ptr<Tensor> rdm0deriv_;
    std::shared_ptr<Tensor> rdm1deriv_;
    std::shared_ptr<Tensor> rdm2deriv_;
    std::shared_ptr<Tensor> rdm3deriv_;
    std::shared_ptr<Tensor> rdm4deriv_;

    std::shared_ptr<Tensor> sigma_;

    // the diagonal denominator
    std::vector<double> eig_;

    // printing functions called from the solve function of a derived class
    static void print_iteration();
    static void print_iteration(const int i, const double en, const double err, const double tim, const int istate = -1);
    static void print_iteration(const bool noconv);

    // compute e0 which is defined as Trace(f(x,x), gamma(x,x))
    double compute_e0();

    // denominator objects
    std::vector<std::shared_ptr<const Denom>> denom_;

    // update t from the residual and denominator (this function does not zero out).
    void update_amplitude(std::shared_ptr<Tensor> t, std::shared_ptr<const Tensor> r, std::shared_ptr<const Denom> denom = nullptr) const;
    void update_amplitude(std::shared_ptr<MultiTensor> t, std::shared_ptr<const MultiTensor> r) const;

    // initialize t2 amplitude
    std::shared_ptr<Tensor> init_amplitude() const;

    // diagonal part of CASPT2 (for efficiency)
    void diagonal(std::shared_ptr<Tensor> r, std::shared_ptr<const Tensor> t) const;

  public:
    SpinFreeMethod(std::shared_ptr<const SMITH_Info> r);

    IndexRange& virt() { return virt_; }
    IndexRange& all() { return all_; }
    IndexRange& closed() { return closed_; }

    std::shared_ptr<const SMITH_Info> ref() const { return ref_; }

    std::shared_ptr<const Civec> civec() const { return civec_; }

    std::shared_ptr<const Coeff> coeff() const { return coeff_; }

    double e0() const { return e0_; }

    double energy() const { return energy_; }

    virtual void solve() = 0;

    std::shared_ptr<const Civec> rdm0deriv() const {
      return rdm0deriv_->civec(det_);
    }

    double dot_product_transpose(std::shared_ptr<const Tensor> r, std::shared_ptr<const Tensor> t2) const;
    double dot_product_transpose(std::shared_ptr<const MultiTensor> r, std::shared_ptr<const MultiTensor> t2) const;
};

}
}

#endif
