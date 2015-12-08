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
    std::shared_ptr<TATensor<std::complex<double>,4>> t2;
    std::shared_ptr<TATensor<std::complex<double>,4>> r;
    std::shared_ptr<TATensor<std::complex<double>,4>> s;
    std::shared_ptr<TATensor<std::complex<double>,4>> n;

    int nstates_;
    std::vector<double> energy_;

    std::vector<std::shared_ptr<MultiTATensor<std::complex<double>,4>>> t2all_;
    std::vector<std::shared_ptr<MultiTATensor<std::complex<double>,4>>> rall_;
    std::vector<std::shared_ptr<MultiTATensor<std::complex<double>,4>>> sall_;
    std::vector<std::shared_ptr<MultiTATensor<std::complex<double>,4>>> nall_;
    void diagonal(std::shared_ptr<TATensor<std::complex<double>,4>> r, std::shared_ptr<const TATensor<std::complex<double>,4>> t) const;


    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma0_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma1_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma2_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma58_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma59_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma60_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma63_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma64_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma65_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma66_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma67_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma3_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma4_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma5_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma74_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma75_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma77_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma78_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma79_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma81_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma84_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma86_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma92_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma95_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma98_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma411_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma412_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma9_();
    std::shared_ptr<FutureTATensor<std::complex<double>,2>> Gamma11_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma160_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma99_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma105_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma110_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma128_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma151_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma159_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma174_();
    std::shared_ptr<FutureTATensor<std::complex<double>,2>> Gamma413_();
    std::shared_ptr<FutureTATensor<std::complex<double>,2>> Gamma415_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma24_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma26_();
    std::shared_ptr<FutureTATensor<std::complex<double>,2>> Gamma27_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma179_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma183_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma193_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma203_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma204_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma229_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma417_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma418_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma31_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma32_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma33_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma230_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma231_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma232_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma233_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma236_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma237_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma238_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma245_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma246_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma253_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma419_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma420_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma317_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma318_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma335_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma336_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma368_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma363_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma389_();
    std::shared_ptr<FutureTATensor<std::complex<double>,2>> Gamma425_();
    std::shared_ptr<FutureTATensor<std::complex<double>,2>> Gamma427_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma394_();
    std::shared_ptr<FutureTATensor<std::complex<double>,8>> Gamma395_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma407_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma409_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma410_();
    std::shared_ptr<FutureTATensor<std::complex<double>,0>> Gamma421_();
    std::shared_ptr<FutureTATensor<std::complex<double>,0>> Gamma423_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma429_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma430_();
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

