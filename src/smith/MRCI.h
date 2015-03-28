//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI.h
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


#ifndef __SRC_SMITH_MRCI_H
#define __SRC_SMITH_MRCI_H

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
namespace MRCI{

class MRCI : public SpinFreeMethod {
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
    std::shared_ptr<FutureTensor> Gamma80_();
    std::shared_ptr<FutureTensor> Gamma81_();
    std::shared_ptr<FutureTensor> Gamma82_();
    std::shared_ptr<FutureTensor> Gamma85_();
    std::shared_ptr<FutureTensor> Gamma86_();
    std::shared_ptr<FutureTensor> Gamma87_();
    std::shared_ptr<FutureTensor> Gamma88_();
    std::shared_ptr<FutureTensor> Gamma89_();
    std::shared_ptr<FutureTensor> Gamma94_();
    std::shared_ptr<FutureTensor> Gamma3_();
    std::shared_ptr<FutureTensor> Gamma4_();
    std::shared_ptr<FutureTensor> Gamma5_();
    std::shared_ptr<FutureTensor> Gamma7_();
    std::shared_ptr<FutureTensor> Gamma97_();
    std::shared_ptr<FutureTensor> Gamma98_();
    std::shared_ptr<FutureTensor> Gamma100_();
    std::shared_ptr<FutureTensor> Gamma101_();
    std::shared_ptr<FutureTensor> Gamma102_();
    std::shared_ptr<FutureTensor> Gamma104_();
    std::shared_ptr<FutureTensor> Gamma107_();
    std::shared_ptr<FutureTensor> Gamma109_();
    std::shared_ptr<FutureTensor> Gamma114_();
    std::shared_ptr<FutureTensor> Gamma115_();
    std::shared_ptr<FutureTensor> Gamma119_();
    std::shared_ptr<FutureTensor> Gamma122_();
    std::shared_ptr<FutureTensor> Gamma550_();
    std::shared_ptr<FutureTensor> Gamma551_();
    std::shared_ptr<FutureTensor> Gamma10_();
    std::shared_ptr<FutureTensor> Gamma12_();
    std::shared_ptr<FutureTensor> Gamma18_();
    std::shared_ptr<FutureTensor> Gamma197_();
    std::shared_ptr<FutureTensor> Gamma126_();
    std::shared_ptr<FutureTensor> Gamma132_();
    std::shared_ptr<FutureTensor> Gamma137_();
    std::shared_ptr<FutureTensor> Gamma155_();
    std::shared_ptr<FutureTensor> Gamma176_();
    std::shared_ptr<FutureTensor> Gamma178_();
    std::shared_ptr<FutureTensor> Gamma179_();
    std::shared_ptr<FutureTensor> Gamma196_();
    std::shared_ptr<FutureTensor> Gamma552_();
    std::shared_ptr<FutureTensor> Gamma554_();
    std::shared_ptr<FutureTensor> Gamma24_();
    std::shared_ptr<FutureTensor> Gamma25_();
    std::shared_ptr<FutureTensor> Gamma27_();
    std::shared_ptr<FutureTensor> Gamma29_();
    std::shared_ptr<FutureTensor> Gamma31_();
    std::shared_ptr<FutureTensor> Gamma32_();
    std::shared_ptr<FutureTensor> Gamma215_();
    std::shared_ptr<FutureTensor> Gamma216_();
    std::shared_ptr<FutureTensor> Gamma217_();
    std::shared_ptr<FutureTensor> Gamma220_();
    std::shared_ptr<FutureTensor> Gamma222_();
    std::shared_ptr<FutureTensor> Gamma221_();
    std::shared_ptr<FutureTensor> Gamma230_();
    std::shared_ptr<FutureTensor> Gamma232_();
    std::shared_ptr<FutureTensor> Gamma234_();
    std::shared_ptr<FutureTensor> Gamma233_();
    std::shared_ptr<FutureTensor> Gamma235_();
    std::shared_ptr<FutureTensor> Gamma240_();
    std::shared_ptr<FutureTensor> Gamma244_();
    std::shared_ptr<FutureTensor> Gamma250_();
    std::shared_ptr<FutureTensor> Gamma251_();
    std::shared_ptr<FutureTensor> Gamma252_();
    std::shared_ptr<FutureTensor> Gamma276_();
    std::shared_ptr<FutureTensor> Gamma568_();
    std::shared_ptr<FutureTensor> Gamma569_();
    std::shared_ptr<FutureTensor> Gamma572_();
    std::shared_ptr<FutureTensor> Gamma573_();
    std::shared_ptr<FutureTensor> Gamma278_();
    std::shared_ptr<FutureTensor> Gamma296_();
    std::shared_ptr<FutureTensor> Gamma312_();
    std::shared_ptr<FutureTensor> Gamma313_();
    std::shared_ptr<FutureTensor> Gamma338_();
    std::shared_ptr<FutureTensor> Gamma48_();
    std::shared_ptr<FutureTensor> Gamma49_();
    std::shared_ptr<FutureTensor> Gamma50_();
    std::shared_ptr<FutureTensor> Gamma51_();
    std::shared_ptr<FutureTensor> Gamma339_();
    std::shared_ptr<FutureTensor> Gamma340_();
    std::shared_ptr<FutureTensor> Gamma341_();
    std::shared_ptr<FutureTensor> Gamma342_();
    std::shared_ptr<FutureTensor> Gamma345_();
    std::shared_ptr<FutureTensor> Gamma346_();
    std::shared_ptr<FutureTensor> Gamma349_();
    std::shared_ptr<FutureTensor> Gamma350_();
    std::shared_ptr<FutureTensor> Gamma351_();
    std::shared_ptr<FutureTensor> Gamma359_();
    std::shared_ptr<FutureTensor> Gamma366_();
    std::shared_ptr<FutureTensor> Gamma556_();
    std::shared_ptr<FutureTensor> Gamma557_();
    std::shared_ptr<FutureTensor> Gamma471_();
    std::shared_ptr<FutureTensor> Gamma503_();
    std::shared_ptr<FutureTensor> Gamma526_();
    std::shared_ptr<FutureTensor> Gamma562_();
    std::shared_ptr<FutureTensor> Gamma564_();
    std::shared_ptr<FutureTensor> Gamma531_();
    std::shared_ptr<FutureTensor> Gamma532_();
    std::shared_ptr<FutureTensor> Gamma533_();
    std::shared_ptr<FutureTensor> Gamma545_();
    std::shared_ptr<FutureTensor> Gamma548_();
    std::shared_ptr<FutureTensor> Gamma549_();
    std::shared_ptr<FutureTensor> Gamma558_();
    std::shared_ptr<FutureTensor> Gamma560_();
    std::shared_ptr<FutureTensor> Gamma566_();
    std::shared_ptr<FutureTensor> Gamma567_();
    std::shared_ptr<Queue> make_residualq(const bool reset = true);
    std::shared_ptr<Queue> make_sourceq(const bool reset = true);
    std::shared_ptr<Queue> make_normq(const bool reset = true);

  public:
    MRCI(std::shared_ptr<const SMITH_Info> ref);
    ~MRCI() {}

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

