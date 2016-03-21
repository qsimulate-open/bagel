//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: caspt2grad.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __SRC_SMITH_CASPT2GRAD_H
#define __SRC_SMITH_CASPT2GRAD_H

#include <src/wfn/reference.h>
#include <src/ci/fci/fci.h>
#include <src/smith/tensor.h>

namespace bagel {

class CASPT2Grad : public Method {
  public:
    using Tensor = SMITH::Tensor_<double>;
  protected:
    std::shared_ptr<const Matrix> coeff_;
    // second-order density matrix
    std::shared_ptr<const Matrix> d1_;
    // first-order density matrix
    std::shared_ptr<const Matrix> d11_;
    // two-body first-order density matrix
    std::shared_ptr<const Tensor> d2_;

    // y from SMITH code
    std::shared_ptr<Civec> cideriv_;
    // FCI utility
    std::shared_ptr<FCI> fci_;

    // for gradient
    int target_;
    int ncore_;
    double energy_;
    double thresh_;

    // properties
    bool do_hyperfine_;

    std::vector<double> ref_energy_;

    std::shared_ptr<DFFullDist> contract_D1(std::shared_ptr<const DFFullDist> full) const;

  public:
    CASPT2Grad(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    void compute() override;

    std::shared_ptr<const Matrix> coeff() const { return coeff_; }
    std::shared_ptr<const Matrix> d1() const { return d1_; }
    std::shared_ptr<const Matrix> d11() const { return d11_; }
    std::shared_ptr<const Tensor> d2() const { return d2_; }

    std::shared_ptr<const Civec> cideriv() const { return cideriv_; }
    std::shared_ptr<FCI> fci() const { return fci_; }

    int target() const { return target_; }
    int ncore() const { return ncore_; }
    double energy() const { return energy_; }
    double thresh() const { return thresh_; }

    bool do_hyperfine() const { return do_hyperfine_; }

    std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; }

    std::tuple<std::shared_ptr<Matrix>,std::shared_ptr<const DFFullDist>>
      compute_Y(std::shared_ptr<const DFHalfDist> half, std::shared_ptr<const DFHalfDist> halfj, std::shared_ptr<const DFHalfDist> halfjj);

    std::shared_ptr<Matrix> diagonal_D1() const;
    std::shared_ptr<Matrix> spin_density_unrelaxed() const;
    std::shared_ptr<Matrix> spin_density_relax(std::shared_ptr<const RDM<1>> zrdm1, std::shared_ptr<const RDM<2>> zrdm2, std::shared_ptr<const Matrix> zmat) const;
};

}

#endif
