//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: SPCASPT2.h
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


#ifndef __SRC_SMITH_SPCASPT2_H
#define __SRC_SMITH_SPCASPT2_H

#include <src/smith/caspt2/CASPT2.h>

namespace bagel {
namespace SMITH {
namespace SPCASPT2{

class SPCASPT2 {
  protected:
    std::shared_ptr<const SMITH_Info<double>> info_;

    IndexRange virt_;
    IndexRange active_;
    IndexRange closed_;
    std::shared_ptr<const IndexRange> rvirt_;
    std::shared_ptr<const IndexRange> ractive_;
    std::shared_ptr<const IndexRange> rclosed_;

    std::shared_ptr<Tensor> t2;

    std::shared_ptr<Tensor> den1;
    std::shared_ptr<Tensor> den2;

    // contains the current RDMs to be used in smith
    std::shared_ptr<Tensor> rdm0_;
    std::shared_ptr<Tensor> rdm1_;
    std::shared_ptr<Tensor> rdm2_;
    std::shared_ptr<Tensor> rdm3_;
    std::shared_ptr<Tensor> rdm4_;

    // alpha density matrices
    std::shared_ptr<Tensor> ardm0_;
    std::shared_ptr<Tensor> ardm1_;
    std::shared_ptr<Tensor> ardm2_;
    std::shared_ptr<Tensor> ardm3_;
    std::shared_ptr<Tensor> ardm4_;

    std::shared_ptr<FutureTensor> Gamma0_();
    std::shared_ptr<FutureTensor> Gamma31_();
    std::shared_ptr<FutureTensor> Gamma34_();
    std::shared_ptr<FutureTensor> Gamma92_();
    std::shared_ptr<FutureTensor> Gamma1_();
    std::shared_ptr<FutureTensor> Gamma32_();
    std::shared_ptr<FutureTensor> Gamma35_();
    std::shared_ptr<FutureTensor> Gamma2_();
    std::shared_ptr<FutureTensor> Gamma37_();
    std::shared_ptr<FutureTensor> Gamma3_();
    std::shared_ptr<FutureTensor> Gamma40_();
    std::shared_ptr<FutureTensor> Gamma4_();
    std::shared_ptr<FutureTensor> Gamma56_();
    std::shared_ptr<FutureTensor> Gamma57_();
    std::shared_ptr<FutureTensor> Gamma5_();
    std::shared_ptr<FutureTensor> Gamma58_();
    std::shared_ptr<FutureTensor> Gamma6_();
    std::shared_ptr<FutureTensor> Gamma7_();
    std::shared_ptr<FutureTensor> Gamma8_();
    std::shared_ptr<FutureTensor> Gamma60_();
    std::shared_ptr<FutureTensor> Gamma61_();
    std::shared_ptr<FutureTensor> Gamma9_();
    std::shared_ptr<FutureTensor> Gamma62_();
    std::shared_ptr<FutureTensor> Gamma11_();
    std::shared_ptr<FutureTensor> Gamma12_();
    std::shared_ptr<FutureTensor> Gamma13_();
    std::shared_ptr<FutureTensor> Gamma65_();
    std::shared_ptr<FutureTensor> Gamma67_();
    std::shared_ptr<FutureTensor> Gamma14_();
    std::shared_ptr<FutureTensor> Gamma81_();
    std::shared_ptr<FutureTensor> Gamma16_();
    std::shared_ptr<FutureTensor> Gamma17_();
    std::shared_ptr<FutureTensor> Gamma22_();
    std::shared_ptr<FutureTensor> Gamma28_();
    std::shared_ptr<FutureTensor> Gamma29_();
    std::shared_ptr<FutureTensor> Gamma33_();
    std::shared_ptr<FutureTensor> Gamma90_();
    std::shared_ptr<FutureTensor> Gamma49_();
    std::shared_ptr<FutureTensor> Gamma51_();
    std::shared_ptr<FutureTensor> Gamma59_();
    std::shared_ptr<Queue> make_densityq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_density1q(const bool reset = true, const bool diagonal = true);

  public:
    SPCASPT2(const CASPT2::CASPT2&);
    ~SPCASPT2() {}

    void solve();

    std::shared_ptr<const Matrix> rdm11() const { return den1->matrix(); }
    std::shared_ptr<const Matrix> rdm12() const { return den2->matrix(); }
};

}
}
}
#endif

