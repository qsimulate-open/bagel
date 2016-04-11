//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2.h
// Copyright (C) 2014 Shiozaki group
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
#include <src/scf/hf/fock.h>
#include <src/util/f77.h>
#include <src/smith/queue.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {
namespace SPCASPT2 { class SPCASPT2; }

namespace CASPT2{

class CASPT2 : public SpinFreeMethod<double> {
  friend class SPCASPT2::SPCASPT2;
  protected:
    std::shared_ptr<Tensor> t2;
    std::shared_ptr<Tensor> r;
    std::shared_ptr<Tensor> s;
    std::shared_ptr<Tensor> n;

    int nstates_;
    std::vector<double> err_;
    std::vector<double> pt2energy_;
    std::vector<std::shared_ptr<MultiTensor>> t2all_;
    std::vector<std::shared_ptr<MultiTensor>> rall_;
    std::vector<std::shared_ptr<MultiTensor>> sall_;
    std::vector<std::shared_ptr<MultiTensor>> lall_;

    std::shared_ptr<Tensor> den1;
    std::shared_ptr<Tensor> den2;
    std::shared_ptr<Tensor> Den1;
    double correlated_norm_;
    std::shared_ptr<Tensor> deci;
    std::shared_ptr<Civec> ci_deriv_;

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
    std::shared_ptr<FutureTensor> Gamma248_();
    std::shared_ptr<FutureTensor> Gamma249_();
    std::shared_ptr<FutureTensor> Gamma250_();
    std::shared_ptr<FutureTensor> Gamma251_();
    std::shared_ptr<FutureTensor> Gamma252_();
    std::shared_ptr<FutureTensor> Gamma253_();
    std::shared_ptr<FutureTensor> Gamma254_();
    std::shared_ptr<FutureTensor> Gamma255_();
    std::shared_ptr<FutureTensor> Gamma257_();
    std::shared_ptr<FutureTensor> Gamma260_();
    std::shared_ptr<FutureTensor> Gamma262_();
    std::shared_ptr<FutureTensor> Gamma264_();
    std::shared_ptr<FutureTensor> Gamma270_();
    std::shared_ptr<FutureTensor> Gamma276_();
    std::shared_ptr<FutureTensor> Gamma277_();
    std::shared_ptr<FutureTensor> Gamma279_();
    std::shared_ptr<FutureTensor> Gamma280_();
    std::shared_ptr<FutureTensor> Gamma282_();
    std::shared_ptr<FutureTensor> Gamma283_();
    std::shared_ptr<FutureTensor> Gamma285_();
    std::shared_ptr<FutureTensor> Gamma286_();
    std::shared_ptr<FutureTensor> Gamma299_();
    std::shared_ptr<FutureTensor> Gamma304_();
    std::shared_ptr<FutureTensor> Gamma305_();
    std::shared_ptr<FutureTensor> Gamma306_();
    std::shared_ptr<FutureTensor> Gamma307_();
    std::shared_ptr<FutureTensor> Gamma308_();
    std::shared_ptr<FutureTensor> Gamma317_();
    std::shared_ptr<FutureTensor> Gamma329_();
    std::shared_ptr<FutureTensor> Gamma340_();
    std::shared_ptr<FutureTensor> Gamma355_();
    std::shared_ptr<FutureTensor> Gamma373_();

    std::shared_ptr<Queue> make_residualq(const bool reset = true, const bool diagonal = true);
    void make_residualq1(std::shared_ptr<Queue>, std::shared_ptr<Task>, const bool);
    void make_residualq2(std::shared_ptr<Queue>, std::shared_ptr<Task>, const bool);
    void make_residualq3(std::shared_ptr<Queue>, std::shared_ptr<Task>, const bool);

    std::shared_ptr<Queue> make_sourceq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_normq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_density1q(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_density2q(const bool reset = true, const bool diagonal = true);

    std::shared_ptr<Queue> make_densityq(const bool reset = true, const bool diagonal = true);
    void make_densityq1(std::shared_ptr<Queue>, std::shared_ptr<Task>, const bool);
    void make_densityq2(std::shared_ptr<Queue>, std::shared_ptr<Task>, const bool);
    void make_densityq3(std::shared_ptr<Queue>, std::shared_ptr<Task>, const bool);

    std::shared_ptr<Queue> make_deciq(const bool reset = true, const bool diagonal = true);
    void make_deciq1(std::shared_ptr<Queue>, std::shared_ptr<Task>, std::shared_ptr<Task>, const bool, std::shared_ptr<Tensor>);
    void make_deciq2(std::shared_ptr<Queue>, std::shared_ptr<Task>, std::shared_ptr<Task>, const bool, std::shared_ptr<Tensor>);
    void make_deciq3(std::shared_ptr<Queue>, std::shared_ptr<Task>, std::shared_ptr<Task>, const bool, std::shared_ptr<Tensor>);

    std::vector<std::shared_ptr<MultiTensor_<double>>>
      solve_linear(std::vector<std::shared_ptr<MultiTensor_<double>>> s, std::vector<std::shared_ptr<MultiTensor_<double>>> t);

  public:
    CASPT2(std::shared_ptr<const SMITH_Info<double>> ref);
    ~CASPT2() {}

    void solve();
    void solve_deriv();

    std::shared_ptr<const Matrix> rdm11() const { return den1->matrix(); }
    std::shared_ptr<const Matrix> rdm12() const { return den2->matrix(); }
    std::shared_ptr<const Tensor> rdm21() const { return Den1; }

    double correlated_norm() const { return correlated_norm_; }

    std::shared_ptr<const Civec> ci_deriv() const { return ci_deriv_; }

};

}
}
}
#endif

