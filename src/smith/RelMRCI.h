//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI.h
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


#ifndef __SRC_SMITH_RelMRCI_H
#define __SRC_SMITH_RelMRCI_H

#include <iostream>
#include <tuple>
#include <iomanip>
#include <src/smith/spinfreebase.h>
#include <src/smith/futuretensor.h>
#include <src/scf/hf/fock.h>
#include <src/util/f77.h>
#include <src/smith/queue.h>
#include <src/smith/multitensor.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class RelMRCI : public SpinFreeMethod<std::complex<double>> {
  protected:
    std::shared_ptr<Tensor> t2;
    std::shared_ptr<Tensor> r;
    std::shared_ptr<Tensor> s;
    std::shared_ptr<Tensor> n;

    int nstates_;
    std::vector<double> energy_;

    std::vector<std::shared_ptr<MultiTensor>> t2all_;
    std::vector<std::shared_ptr<MultiTensor>> rall_;
    std::vector<std::shared_ptr<MultiTensor>> sall_;
    std::vector<std::shared_ptr<MultiTensor>> nall_;


    std::shared_ptr<FutureTensor> Gamma0_();
    std::shared_ptr<FutureTensor> Gamma1_();
    std::shared_ptr<FutureTensor> Gamma2_();
    std::shared_ptr<FutureTensor> Gamma84_();
    std::shared_ptr<FutureTensor> Gamma85_();
    std::shared_ptr<FutureTensor> Gamma86_();
    std::shared_ptr<FutureTensor> Gamma89_();
    std::shared_ptr<FutureTensor> Gamma90_();
    std::shared_ptr<FutureTensor> Gamma91_();
    std::shared_ptr<FutureTensor> Gamma92_();
    std::shared_ptr<FutureTensor> Gamma93_();
    std::shared_ptr<FutureTensor> Gamma98_();
    std::shared_ptr<FutureTensor> Gamma3_();
    std::shared_ptr<FutureTensor> Gamma4_();
    std::shared_ptr<FutureTensor> Gamma5_();
    std::shared_ptr<FutureTensor> Gamma7_();
    std::shared_ptr<FutureTensor> Gamma101_();
    std::shared_ptr<FutureTensor> Gamma102_();
    std::shared_ptr<FutureTensor> Gamma104_();
    std::shared_ptr<FutureTensor> Gamma105_();
    std::shared_ptr<FutureTensor> Gamma106_();
    std::shared_ptr<FutureTensor> Gamma108_();
    std::shared_ptr<FutureTensor> Gamma111_();
    std::shared_ptr<FutureTensor> Gamma113_();
    std::shared_ptr<FutureTensor> Gamma118_();
    std::shared_ptr<FutureTensor> Gamma119_();
    std::shared_ptr<FutureTensor> Gamma123_();
    std::shared_ptr<FutureTensor> Gamma126_();
    std::shared_ptr<FutureTensor> Gamma558_();
    std::shared_ptr<FutureTensor> Gamma559_();
    std::shared_ptr<FutureTensor> Gamma10_();
    std::shared_ptr<FutureTensor> Gamma12_();
    std::shared_ptr<FutureTensor> Gamma18_();
    std::shared_ptr<FutureTensor> Gamma201_();
    std::shared_ptr<FutureTensor> Gamma130_();
    std::shared_ptr<FutureTensor> Gamma136_();
    std::shared_ptr<FutureTensor> Gamma141_();
    std::shared_ptr<FutureTensor> Gamma159_();
    std::shared_ptr<FutureTensor> Gamma180_();
    std::shared_ptr<FutureTensor> Gamma182_();
    std::shared_ptr<FutureTensor> Gamma183_();
    std::shared_ptr<FutureTensor> Gamma200_();
    std::shared_ptr<FutureTensor> Gamma560_();
    std::shared_ptr<FutureTensor> Gamma562_();
    std::shared_ptr<FutureTensor> Gamma24_();
    std::shared_ptr<FutureTensor> Gamma25_();
    std::shared_ptr<FutureTensor> Gamma27_();
    std::shared_ptr<FutureTensor> Gamma28_();
    std::shared_ptr<FutureTensor> Gamma30_();
    std::shared_ptr<FutureTensor> Gamma31_();
    std::shared_ptr<FutureTensor> Gamma33_();
    std::shared_ptr<FutureTensor> Gamma34_();
    std::shared_ptr<FutureTensor> Gamma219_();
    std::shared_ptr<FutureTensor> Gamma220_();
    std::shared_ptr<FutureTensor> Gamma221_();
    std::shared_ptr<FutureTensor> Gamma224_();
    std::shared_ptr<FutureTensor> Gamma226_();
    std::shared_ptr<FutureTensor> Gamma225_();
    std::shared_ptr<FutureTensor> Gamma234_();
    std::shared_ptr<FutureTensor> Gamma235_();
    std::shared_ptr<FutureTensor> Gamma237_();
    std::shared_ptr<FutureTensor> Gamma239_();
    std::shared_ptr<FutureTensor> Gamma238_();
    std::shared_ptr<FutureTensor> Gamma240_();
    std::shared_ptr<FutureTensor> Gamma245_();
    std::shared_ptr<FutureTensor> Gamma246_();
    std::shared_ptr<FutureTensor> Gamma250_();
    std::shared_ptr<FutureTensor> Gamma256_();
    std::shared_ptr<FutureTensor> Gamma257_();
    std::shared_ptr<FutureTensor> Gamma258_();
    std::shared_ptr<FutureTensor> Gamma282_();
    std::shared_ptr<FutureTensor> Gamma284_();
    std::shared_ptr<FutureTensor> Gamma303_();
    std::shared_ptr<FutureTensor> Gamma320_();
    std::shared_ptr<FutureTensor> Gamma321_();
    std::shared_ptr<FutureTensor> Gamma346_();
    std::shared_ptr<FutureTensor> Gamma52_();
    std::shared_ptr<FutureTensor> Gamma53_();
    std::shared_ptr<FutureTensor> Gamma54_();
    std::shared_ptr<FutureTensor> Gamma55_();
    std::shared_ptr<FutureTensor> Gamma347_();
    std::shared_ptr<FutureTensor> Gamma348_();
    std::shared_ptr<FutureTensor> Gamma349_();
    std::shared_ptr<FutureTensor> Gamma350_();
    std::shared_ptr<FutureTensor> Gamma353_();
    std::shared_ptr<FutureTensor> Gamma354_();
    std::shared_ptr<FutureTensor> Gamma357_();
    std::shared_ptr<FutureTensor> Gamma358_();
    std::shared_ptr<FutureTensor> Gamma359_();
    std::shared_ptr<FutureTensor> Gamma367_();
    std::shared_ptr<FutureTensor> Gamma374_();
    std::shared_ptr<FutureTensor> Gamma564_();
    std::shared_ptr<FutureTensor> Gamma565_();
    std::shared_ptr<FutureTensor> Gamma479_();
    std::shared_ptr<FutureTensor> Gamma511_();
    std::shared_ptr<FutureTensor> Gamma534_();
    std::shared_ptr<FutureTensor> Gamma570_();
    std::shared_ptr<FutureTensor> Gamma572_();
    std::shared_ptr<FutureTensor> Gamma539_();
    std::shared_ptr<FutureTensor> Gamma540_();
    std::shared_ptr<FutureTensor> Gamma541_();
    std::shared_ptr<FutureTensor> Gamma553_();
    std::shared_ptr<FutureTensor> Gamma556_();
    std::shared_ptr<FutureTensor> Gamma557_();
    std::shared_ptr<FutureTensor> Gamma566_();
    std::shared_ptr<FutureTensor> Gamma568_();
    std::shared_ptr<FutureTensor> Gamma574_();
    std::shared_ptr<FutureTensor> Gamma575_();
    std::shared_ptr<Queue> make_residualq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_sourceq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_normq(const bool reset = true, const bool diagonal = true);

  public:
    RelMRCI(std::shared_ptr<const SMITH_Info<std::complex<double>>> ref);
    ~RelMRCI() {}

    void solve();
    void solve_deriv();

    double accumulate(std::shared_ptr<Queue> queue) {
      double sum = 0.0;
      while (!queue->done())
        sum += queue->next_compute()->target();
      return sum;
    }

};

}
}
}
#endif

