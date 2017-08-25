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
#include <src/ci/fci/distfci.h>
#include <src/smith/tensor.h>
#include <src/smith/smith.h>
#include <src/grad/nacmtype.h>

namespace bagel {

class CASPT2Grad : public Method {
  public:
    using Tensor = SMITH::Tensor_<double>;

  protected:
    // second-order density matrix
    std::shared_ptr<const Matrix> d1_;
    std::shared_ptr<Matrix> vd1_;
    // XMS density if available
    std::shared_ptr<const Matrix> dcheck_;
    double energy_;
    std::vector<double> dipole_;

    std::shared_ptr<const Matrix> coeff_;
    // first-order density matrix
    std::shared_ptr<const Matrix> d11_;
    std::shared_ptr<const Tensor> d2_;
    // zeroth-order density matrix (mixed state)
    std::shared_ptr<RDM<1>> d10ms_;
    std::shared_ptr<RDM<2>> d20ms_;
    // norm of the first-order wave function
    std::vector<double> wf1norm_;

    // second-order spin density matrix
    std::shared_ptr<const Matrix> sd1_;
    // first-order spin density matrix
    std::shared_ptr<const Matrix> sd11_;
    // rotation matrix in MS-CASPT2
    std::shared_ptr<const Matrix> msrot_;

    // y from SMITH code
    std::shared_ptr<Dvec> cideriv_;
    // FCI utility
    std::shared_ptr<FCI_base> fci_;

    // SMITH
    std::shared_ptr<Smith> smith_;

    // for gradient
    int nstates_;
    int ncore_;
    int target_;
    double thresh_;

    int maxziter_;

    // properties
    bool do_hyperfine_;
    std::shared_ptr<DFFullDist> contract_D1(std::shared_ptr<const DFFullDist> full) const;

  public:
    CASPT2Grad(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    const double& msrot(int i, int j) const { return msrot_->element(i, j); }
    std::shared_ptr<const Matrix> coeff() const { return coeff_; }
    std::shared_ptr<const Matrix> d11() const { return d11_; }
    std::shared_ptr<const Tensor> d2() const { return d2_; }
    std::shared_ptr<const RDM<1>> d10ms() const { return d10ms_; }
    std::shared_ptr<const RDM<2>> d20ms() const { return d20ms_; }

    std::shared_ptr<const Dvec> cideriv() const { return cideriv_; }
    std::shared_ptr<FCI_base> fci() const { return fci_; }
    std::shared_ptr<Smith> smith() const { return smith_; }

    int nstates() const { return nstates_; }
    int ncore() const { return ncore_; }
    int target() const { return target_; }
    double thresh() const { return thresh_; }

    int maxziter() const { return maxziter_; }

    bool do_hyperfine() const { return do_hyperfine_; }

    std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; }

    void compute() override;
    void compute_gradient(const int istate, const int jstate, std::shared_ptr<const NacmType> nacmtype = std::make_shared<const NacmType>("interstate"), const bool nocider = false);

    std::shared_ptr<const Matrix> d1() const { return d1_; }
    std::shared_ptr<const Matrix> vd1() const { return vd1_; }
    std::shared_ptr<const Matrix> dcheck() const { return dcheck_; }
    double energy() const { return energy_; }
    const std::vector<double>& dipole() const { return dipole_; }
    double dipole(const int i) const { return dipole_[i]; }

    // internal functions

    std::shared_ptr<Matrix> diagonal_D1() const;
    std::shared_ptr<Matrix> spin_density_unrelaxed() const;
    std::shared_ptr<Matrix> spin_density_relax(std::shared_ptr<const RDM<1>> zrdm1, std::shared_ptr<const RDM<2>> zrdm2, std::shared_ptr<const Matrix> zmat) const;

    // function for gradient
    std::tuple<std::shared_ptr<Matrix>,std::shared_ptr<const DFFullDist>>
      compute_Y(std::shared_ptr<const DFHalfDist> half, std::shared_ptr<const DFHalfDist> halfj, std::shared_ptr<const DFHalfDist> halfjj, const bool nacme);

    // function for nacme
    void augment_Y(std::shared_ptr<Matrix> d0ms, std::shared_ptr<Matrix> g0, std::shared_ptr<Dvec> g1, std::shared_ptr<const DFHalfDist> halfj, const int istate, const int jstate, const double egap);
};

// Single point calc. (only for finite difference nacme)

class CASPT2Energy : public Method {
  public:
    using Tensor = SMITH::Tensor_<double>;
  protected:
    std::shared_ptr<const Matrix> coeff_;
    std::shared_ptr<const Matrix> msrot_;
    std::vector<double> energy_;
    std::shared_ptr<FCI_base> fci_;
    std::shared_ptr<Matrix> vd1_;
    double thresh_;
    int target_state1_;
    int target_state2_;
    int ncore_;

    bool do_hyperfine_;

  public:
    CASPT2Energy(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    void compute() override;

    const double& msrot(int i, int j) const { return msrot_->element(i, j); }
    std::shared_ptr<const Matrix> coeff() const { return coeff_; }
    std::shared_ptr<const Matrix> msrot() const { return msrot_; }
    std::shared_ptr<const Matrix> vd1() const { return vd1_; }
    int target1() const { return target_state1_; }
    int target2() const { return target_state2_; }
    int ncore() const { return ncore_; }
    double energy(const int i) const { return energy_[i]; }
    const std::vector<double>& energy() const { return energy_; }
    virtual std::shared_ptr<const Reference> conv_to_ref() const override;

};


}

#endif
