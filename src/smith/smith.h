//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: smith.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Matthew K. MacLeod <matthew.macleod@northwestern.edu>
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


#ifndef __SRC_SMITH_SMITH_H
#define __SRC_SMITH_SMITH_H

#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/spinfreebase.h>
#endif
#include <stddef.h>
#include <map>
#include <memory>
#include <src/wfn/method.h>
#include <src/wfn/reference.h>
#include <src/smith/tensor.h>
#include <src/grad/nacmtype.h>

namespace bagel {

class Smith : public Method {
  public:
    using Tensor = SMITH::Tensor_<double>;

  protected:
#ifdef COMPILE_SMITH
    std::shared_ptr<SMITH::SpinFreeMethod<double>> algo_;
#endif

    // correlated density matrices
    // second order density matrix
    std::shared_ptr<const Matrix> dm1_;
    // first order density matrices
    std::shared_ptr<const Matrix> dm11_;
    std::shared_ptr<const Tensor> dm2_;
    // XMS density matrix
    std::shared_ptr<const Matrix> dcheck_;
    // second order spin density matrix
    std::shared_ptr<const Matrix> sdm1_;
    // first order spin density matrix
    std::shared_ptr<const Matrix> sdm11_;
    // norm of the first-order wave function <1|1>
    std::vector<double> wf1norm_;
    // ci derivative
    std::shared_ptr<const Dvec> cider_;
    // rotation matrix in MS-CASPT2
    std::shared_ptr<const Matrix> msrot_;
    std::shared_ptr<const Matrix> vd1_;

    std::shared_ptr<const Matrix> coeff_;

  public:
    Smith(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    void compute() override;
    // Gradient module to be separated
    void compute_gradient(const int istate, const int jstate, std::shared_ptr<const NacmType> nacmtype = std::make_shared<const NacmType>("interstate"), const bool nocider = false);

    // just return the reference used in SMITH code
    std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; }

    std::shared_ptr<const Matrix> dm1() const { return dm1_; }
    std::shared_ptr<const Matrix> dm11() const { return dm11_; }
    std::shared_ptr<const Tensor> dm2() const { return dm2_; }
    std::shared_ptr<const Matrix> dcheck() const { return dcheck_; }
    std::shared_ptr<const Matrix> sdm1() const { return sdm1_; }
    std::shared_ptr<const Matrix> sdm11() const { return sdm11_; }
    std::vector<double> wf1norm() const { return wf1norm_; }
    std::shared_ptr<const Dvec> cideriv() const { return cider_; }
    std::shared_ptr<const Matrix> msrot() const { return msrot_; }
    std::shared_ptr<const Matrix> vd1() const { return vd1_; }

    std::shared_ptr<const Matrix> coeff() const { return coeff_; }

#ifdef COMPILE_SMITH
    std::shared_ptr<const SMITH::SpinFreeMethod<double>> algo() const { return algo_; }
#endif

};


class RelSmith : public Method {
  protected:
#ifdef COMPILE_SMITH
    std::shared_ptr<SMITH::SpinFreeMethod<std::complex<double>>> algo_;
#endif
    std::shared_ptr<const ZMatrix> coeff_;

  public:
    RelSmith(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    void compute() override {
#ifdef COMPILE_SMITH
      algo_->solve();
#endif
    }

    std::shared_ptr<const Reference> conv_to_ref() const override { return std::shared_ptr<const Reference>(); }
    std::shared_ptr<const ZMatrix> coeff() const { return coeff_; }

#ifdef COMPILE_SMITH
    std::shared_ptr<const SMITH::SpinFreeMethod<std::complex<double>>> algo() const { return algo_; }
#endif
};

}
#endif
