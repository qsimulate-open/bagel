//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2.h
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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
namespace CASPT2{

class CASPT2 : public SpinFreeMethod {
  protected:
    using SpinFreeMethod::ref_;
    using SpinFreeMethod::closed_;
    using SpinFreeMethod::active_;
    using SpinFreeMethod::virt_;
    using SpinFreeMethod::ci_;
    using SpinFreeMethod::rclosed_;
    using SpinFreeMethod::ractive_;
    using SpinFreeMethod::rvirt_;
    using SpinFreeMethod::rci_;
    using SpinFreeMethod::h1_;
    using SpinFreeMethod::f1_;
    using SpinFreeMethod::v2_;
    using SpinFreeMethod::rdm1_;
    using SpinFreeMethod::rdm2_;
    using SpinFreeMethod::rdm3_;
    using SpinFreeMethod::rdm4_;
    using SpinFreeMethod::rdm1deriv_;
    using SpinFreeMethod::rdm2deriv_;
    using SpinFreeMethod::rdm3deriv_;
    using SpinFreeMethod::rdm4deriv_;

    std::shared_ptr<Tensor> t2;
    std::shared_ptr<Tensor> r;
    double e0_;
    std::shared_ptr<Tensor> den1;
    std::shared_ptr<Tensor> den2;
    std::shared_ptr<Tensor> Den1;
    double correlated_norm_;
    std::shared_ptr<Tensor> deci;


    std::shared_ptr<FutureTensor> Gamma0_();
    std::shared_ptr<FutureTensor> Gamma94_();
    std::shared_ptr<FutureTensor> Gamma2_();
    std::shared_ptr<FutureTensor> Gamma3_();
    std::shared_ptr<FutureTensor> Gamma4_();
    std::shared_ptr<FutureTensor> Gamma5_();
    std::shared_ptr<FutureTensor> Gamma6_();
    std::shared_ptr<FutureTensor> Gamma7_();
    std::shared_ptr<FutureTensor> Gamma9_();
    std::shared_ptr<FutureTensor> Gamma107_();
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
    std::shared_ptr<FutureTensor> Gamma69_();
    std::shared_ptr<FutureTensor> Gamma81_();
    std::shared_ptr<FutureTensor> Gamma92_();
    std::shared_ptr<FutureTensor> Gamma162_();
    std::shared_ptr<FutureTensor> Gamma193_();
    std::shared_ptr<FutureTensor> Gamma196_();
    std::shared_ptr<FutureTensor> Gamma254_();
    std::shared_ptr<FutureTensor> Gamma167_();
    std::shared_ptr<FutureTensor> Gamma220_();
    std::shared_ptr<FutureTensor> Gamma176_();
    std::shared_ptr<FutureTensor> Gamma272_();
    std::shared_ptr<FutureTensor> Gamma273_();
    std::shared_ptr<FutureTensor> Gamma274_();
    std::shared_ptr<FutureTensor> Gamma275_();
    std::shared_ptr<FutureTensor> Gamma276_();
    std::shared_ptr<FutureTensor> Gamma277_();
    std::shared_ptr<FutureTensor> Gamma278_();
    std::shared_ptr<FutureTensor> Gamma279_();
    std::shared_ptr<FutureTensor> Gamma281_();
    std::shared_ptr<FutureTensor> Gamma284_();
    std::shared_ptr<FutureTensor> Gamma286_();
    std::shared_ptr<FutureTensor> Gamma288_();
    std::shared_ptr<FutureTensor> Gamma294_();
    std::shared_ptr<FutureTensor> Gamma300_();
    std::shared_ptr<FutureTensor> Gamma301_();
    std::shared_ptr<FutureTensor> Gamma303_();
    std::shared_ptr<FutureTensor> Gamma304_();
    std::shared_ptr<FutureTensor> Gamma306_();
    std::shared_ptr<FutureTensor> Gamma307_();
    std::shared_ptr<FutureTensor> Gamma309_();
    std::shared_ptr<FutureTensor> Gamma310_();
    std::shared_ptr<FutureTensor> Gamma323_();
    std::shared_ptr<FutureTensor> Gamma328_();
    std::shared_ptr<FutureTensor> Gamma329_();
    std::shared_ptr<FutureTensor> Gamma330_();
    std::shared_ptr<FutureTensor> Gamma331_();
    std::shared_ptr<FutureTensor> Gamma332_();
    std::shared_ptr<FutureTensor> Gamma341_();
    std::shared_ptr<FutureTensor> Gamma353_();
    std::shared_ptr<FutureTensor> Gamma364_();
    std::shared_ptr<FutureTensor> Gamma379_();
    std::shared_ptr<FutureTensor> Gamma397_();
    std::shared_ptr<Queue> make_residualq();
    std::shared_ptr<Queue> make_energyq();
    std::shared_ptr<Queue> make_corrq();
    std::shared_ptr<Queue> make_densityq();
    std::shared_ptr<Queue> make_density1q();
    std::shared_ptr<Queue> make_density2q();
    std::shared_ptr<Queue> make_deciq();

  public:
    CASPT2(std::shared_ptr<const SMITH_Info> ref);
    ~CASPT2() {}

    void solve();

    double accumulate(std::shared_ptr<Queue> queue) {
      double sum = 0.0;
      while (!queue->done())
        sum += queue->next_compute()->target();
      return sum;
    }

    std::shared_ptr<const Matrix> rdm11() const { return den1->matrix(); }
    std::shared_ptr<const Matrix> rdm12() const { return den2->matrix(); }
    std::shared_ptr<const Matrix> rdm21() const { return Den1->matrix2(); }

    double correlated_norm() const { return correlated_norm_; }

    std::shared_ptr<const Civec> ci_deriv() const { return deci->civec(this->det_); }

};

}
}
}
#endif

