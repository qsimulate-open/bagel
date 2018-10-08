//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software; you can redistribute it and/or modify
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


#ifndef __SRC_SMITH_CASPT2_H
#define __SRC_SMITH_CASPT2_H

#include <iostream>
#include <tuple>
#include <iomanip>
#include <src/smith/spinfreebase.h>
#include <src/smith/futuretensor.h>
#include <src/smith/orthogonal.h>
#include <src/scf/hf/fock.h>
#include <src/util/f77.h>
#include <src/smith/queue.h>
#include <src/smith/smith_info.h>
#include <src/grad/nacmtype.h>

namespace bagel {
namespace SMITH {

namespace SPCASPT2 { class SPCASPT2; }
namespace MSCASPT2 { class MSCASPT2; }

namespace CASPT2{

class CASPT2 : public SpinFreeMethod<double> {
  friend class SPCASPT2::SPCASPT2;
  friend class MSCASPT2::MSCASPT2;
  protected:
    // these are tensors that are used internally
    std::shared_ptr<Tensor> t2;
    std::shared_ptr<Tensor> r;
    std::shared_ptr<Tensor> s;
    std::shared_ptr<Tensor> n;
    std::shared_ptr<Tensor> den1;
    std::shared_ptr<Tensor> den2;
    std::shared_ptr<Tensor> Den1;
    std::shared_ptr<Tensor> den0ci;
    std::shared_ptr<Tensor> den1ci;
    std::shared_ptr<Tensor> den2ci;
    std::shared_ptr<Tensor> den3ci;
    std::shared_ptr<Tensor> den4ci;
    std::shared_ptr<Tensor> den0cit;
    std::shared_ptr<Tensor> den1cit;
    std::shared_ptr<Tensor> den2cit;
    std::shared_ptr<Tensor> den3cit;
    std::shared_ptr<Tensor> den4cit;
    std::shared_ptr<double> den0cirdmt;
    std::shared_ptr<RDM<1>> den1cirdmt;
    std::shared_ptr<RDM<2>> den2cirdmt;
    std::shared_ptr<RDM<3>> den3cirdmt;
    std::shared_ptr<RDM<3>> den4cirdmt;

    int nstates_;
    std::vector<double> err_;
    std::vector<double> pt2energy_;
    std::shared_ptr<Matrix> heff_;

    std::vector<std::shared_ptr<MultiTensor>> t2all_;
    std::vector<std::shared_ptr<MultiTensor>> rall_;
    std::vector<std::shared_ptr<MultiTensor>> sall_;
    std::vector<std::shared_ptr<MultiTensor>> lall_;

    // in orthogonal basis
    std::shared_ptr<Orthogonal_Basis> l_orthogonal_;
    std::shared_ptr<Orthogonal_Basis> t_orthogonal_;

    std::shared_ptr<const Matrix> den1_;
    std::shared_ptr<const Matrix> den2_;
    std::shared_ptr<const Matrix> den2_shift_;
    std::shared_ptr<const Tensor> Den1_;
    // for derivative coupling only
    std::shared_ptr<const Matrix> vden1_;

    std::vector<double> correlated_norm_;
    std::vector<double> correlated_norm_imag_;
    std::shared_ptr<Tensor> deci;
    std::shared_ptr<Dvec> ci_deriv_;
    std::shared_ptr<const Matrix> dcheck_;

    void diagonal(std::shared_ptr<Tensor> r, std::shared_ptr<const Tensor> t, const bool diagonal) const;

    std::shared_ptr<FutureTensor> Gamma0_();
    std::shared_ptr<FutureTensor> Gamma92_();
    std::shared_ptr<FutureTensor> Gamma2_();
    std::shared_ptr<FutureTensor> Gamma3_();
    std::shared_ptr<FutureTensor> Gamma4_();
    std::shared_ptr<FutureTensor> Gamma5_();
    std::shared_ptr<FutureTensor> Gamma6_();
    std::shared_ptr<FutureTensor> Gamma7_();
    std::shared_ptr<FutureTensor> Gamma9_();
    std::shared_ptr<FutureTensor> Gamma12_();
    std::shared_ptr<FutureTensor> Gamma14_();
    std::shared_ptr<FutureTensor> Gamma16_();
    std::shared_ptr<FutureTensor> Gamma22_();
    std::shared_ptr<FutureTensor> Gamma28_();
    std::shared_ptr<FutureTensor> Gamma29_();
    std::shared_ptr<FutureTensor> Gamma31_();
    std::shared_ptr<FutureTensor> Gamma32_();
    std::shared_ptr<FutureTensor> Gamma34_();
    std::shared_ptr<FutureTensor> Gamma35_();
    std::shared_ptr<FutureTensor> Gamma37_();
    std::shared_ptr<FutureTensor> Gamma38_();
    std::shared_ptr<FutureTensor> Gamma51_();
    std::shared_ptr<FutureTensor> Gamma56_();
    std::shared_ptr<FutureTensor> Gamma57_();
    std::shared_ptr<FutureTensor> Gamma58_();
    std::shared_ptr<FutureTensor> Gamma59_();
    std::shared_ptr<FutureTensor> Gamma60_();
    std::shared_ptr<FutureTensor> Gamma79_();
    std::shared_ptr<FutureTensor> Gamma90_();
    std::shared_ptr<FutureTensor> Gamma105_();
    std::shared_ptr<FutureTensor> Gamma138_();
    std::shared_ptr<FutureTensor> Gamma169_();
    std::shared_ptr<FutureTensor> Gamma172_();
    std::shared_ptr<FutureTensor> Gamma230_();
    std::shared_ptr<FutureTensor> Gamma143_();
    std::shared_ptr<FutureTensor> Gamma196_();
    std::shared_ptr<FutureTensor> Gamma152_();

    std::shared_ptr<Queue> make_residualq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_sourceq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_normq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_density1q(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_density2q(const bool reset = true, const bool diagonal = true);

    std::shared_ptr<Queue> make_densityq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_deciq(const bool reset = true, const bool diagonal = true);

    std::vector<std::shared_ptr<MultiTensor_<double>>>
      solve_linear(std::vector<std::shared_ptr<MultiTensor_<double>>> s, std::vector<std::shared_ptr<MultiTensor_<double>>> t);
    std::tuple<std::shared_ptr<Orthogonal_Basis>,std::vector<std::shared_ptr<MultiTensor_<double>>>>
      solve_linear_orthogonal(std::vector<std::shared_ptr<MultiTensor_<double>>> s, std::vector<std::shared_ptr<MultiTensor_<double>>> t);

    std::tuple<std::shared_ptr<double>,std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>,std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<3>>> feed_denci();
    std::shared_ptr<VectorB> contract_rdm_deriv(std::shared_ptr<const CIWfn> ciwfn, int offset, int size, std::shared_ptr<const Matrix> fock);
    void do_rdm_deriv(double factor);

  public:
    CASPT2(std::shared_ptr<const SMITH_Info<double>> ref);
    CASPT2(const CASPT2& cas);
    ~CASPT2() {}

    void solve();
    void solve_gradient(const int targetJ, const int targetI, std::shared_ptr<const NacmType> nacmtype = std::make_shared<const NacmType>(), const bool nocider = false);
    void solve_dm(const int istate, const int jstate);

    std::shared_ptr<const Matrix> msrot() const { return xmsmat_ ? std::make_shared<Matrix>(*xmsmat_ * *heff_) : heff_; }
    std::shared_ptr<const Matrix> rdm11() const { return den1_; }
    std::shared_ptr<const Matrix> rdm12() const { return den2_; }
    std::shared_ptr<const Tensor> rdm21() const { return Den1_; }
    std::shared_ptr<const Matrix> vden1() const { return vden1_; }

    std::vector<double> correlated_norm() const { return correlated_norm_; }

    std::shared_ptr<const Dvec> ci_deriv() const { return ci_deriv_; }
    std::shared_ptr<const Matrix> dcheck() const { return dcheck_; }

};

}
}
}
#endif

