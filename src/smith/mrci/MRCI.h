//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI.h
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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

class MRCI : public SpinFreeMethod<double> {
  protected:
    std::shared_ptr<TATensor<double,4>> t2;
    std::shared_ptr<TATensor<double,4>> r;
    std::shared_ptr<TATensor<double,4>> s;
    std::shared_ptr<TATensor<double,4>> n;

    int nstates_;
    std::vector<double> energy_;

    std::vector<std::shared_ptr<MultiTATensor<double,4>>> t2all_;
    std::vector<std::shared_ptr<MultiTATensor<double,4>>> rall_;
    std::vector<std::shared_ptr<MultiTATensor<double,4>>> sall_;
    std::vector<std::shared_ptr<MultiTATensor<double,4>>> nall_;
    void diagonal(std::shared_ptr<TATensor<double,4>> r, std::shared_ptr<const TATensor<double,4>> t) const;


    std::shared_ptr<FutureTATensor<double,4>> Gamma0_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma1_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma2_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma80_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma81_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma82_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma85_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma86_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma87_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma88_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma89_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma94_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma3_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma4_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma5_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma7_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma97_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma98_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma100_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma101_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma102_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma104_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma107_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma109_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma114_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma115_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma119_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma122_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma547_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma548_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma10_();
    std::shared_ptr<FutureTATensor<double,2>> Gamma12_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma18_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma197_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma126_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma132_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma137_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma155_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma176_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma178_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma179_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma196_();
    std::shared_ptr<FutureTATensor<double,2>> Gamma549_();
    std::shared_ptr<FutureTATensor<double,2>> Gamma551_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma24_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma25_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma27_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma29_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma31_();
    std::shared_ptr<FutureTATensor<double,2>> Gamma32_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma215_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma216_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma217_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma220_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma222_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma221_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma230_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma232_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma234_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma233_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma235_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma240_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma244_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma250_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma251_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma252_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma276_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma565_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma566_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma569_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma570_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma278_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma296_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma312_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma313_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma338_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma48_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma49_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma50_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma51_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma339_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma340_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma341_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma342_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma345_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma346_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma349_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma350_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma351_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma359_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma366_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma553_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma554_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma471_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma503_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma524_();
    std::shared_ptr<FutureTATensor<double,2>> Gamma559_();
    std::shared_ptr<FutureTATensor<double,2>> Gamma561_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma529_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma530_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma531_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma543_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma545_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma546_();
    std::shared_ptr<FutureTATensor<double,0>> Gamma555_();
    std::shared_ptr<FutureTATensor<double,0>> Gamma557_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma563_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma564_();
    std::shared_ptr<Queue> make_residualq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_sourceq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_normq(const bool reset = true, const bool diagonal = true);

  public:
    MRCI(std::shared_ptr<const SMITH_Info<double>> ref);
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

