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
#include <src/smith/denom.h>
#include <src/smith/tensor.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {

class SpinFreeMethod {
  protected:
    // deprecated
    IndexRange virt_;
    IndexRange active_;
    IndexRange closed_;
    IndexRange all_;
    IndexRange ci_;

    // new
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
    std::shared_ptr<Tensor> rdm1_;
    std::shared_ptr<Tensor> rdm2_;
    std::shared_ptr<Tensor> rdm3_;
    std::shared_ptr<Tensor> rdm4_;

    // original determinants (for use in output)
    std::shared_ptr<const Determinants> det_;

    // rdm ci derivatives
    std::shared_ptr<Tensor> rdm0deriv_;
    std::shared_ptr<Tensor> rdm1deriv_;
    std::shared_ptr<Tensor> rdm2deriv_;
    std::shared_ptr<Tensor> rdm3deriv_;
    std::shared_ptr<Tensor> rdm4deriv_;

    std::shared_ptr<Tensor> sigma_;

    mutable std::chrono::high_resolution_clock::time_point time_;

    // the diagonal denominator
    std::vector<double> eig_;

    // printing functions called from the solve function of a derived class
    void print_iteration() const;
    void print_iteration(const int i, const double en, const double err) const;
    void print_iteration(const bool noconv) const;


    // E0 is defined as Trace(f(x,x), gamma(x,x))
    // For instance, E0 is 0 for MP2.
    double compute_e0() const;

    std::shared_ptr<const Denom> denom_;
    std::shared_ptr<const Matrix> shalf_xhh() const { return denom_->shalf_xhh(); }
    std::shared_ptr<const Matrix> shalf_xxh() const { return denom_->shalf_xxh(); }
    std::shared_ptr<const Matrix> shalf_xh() const { return denom_->shalf_xh(); }
    std::shared_ptr<const Matrix> shalf_hh() const { return denom_->shalf_hh(); }
    std::shared_ptr<const Matrix> shalf_xx() const { return denom_->shalf_xx(); }
    std::shared_ptr<const Matrix> shalf_h() const { return denom_->shalf_h(); }
    std::shared_ptr<const Matrix> shalf_x() const { return denom_->shalf_x(); }
    const double& denom_xhh(const size_t i) const { return denom_->denom_xhh(i); }
    const double& denom_xxh(const size_t i) const { return denom_->denom_xxh(i); }
    const double& denom_xh(const size_t i) const { return denom_->denom_xh(i); }
    const double& denom_hh(const size_t i) const { return denom_->denom_hh(i); }
    const double& denom_xx(const size_t i) const { return denom_->denom_xx(i); }
    const double& denom_h(const size_t i) const { return denom_->denom_h(i); }
    const double& denom_x(const size_t i) const { return denom_->denom_x(i); }

    std::shared_ptr<Tensor> init_amplitude() const;
    double dot_product_transpose(std::shared_ptr<const Tensor> r, std::shared_ptr<const Tensor> t2) const;
    void update_amplitude(std::shared_ptr<Tensor> t, std::shared_ptr<const Tensor> r) const;
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
};

}
}

#endif
