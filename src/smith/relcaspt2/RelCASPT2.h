//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2.h
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


#ifndef __SRC_SMITH_RelCASPT2_H
#define __SRC_SMITH_RelCASPT2_H

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
namespace RelCASPT2{

class RelCASPT2 : public SpinFreeMethod<std::complex<double>> {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> t2;
    std::shared_ptr<TATensor<std::complex<double>,4>> r;
    std::shared_ptr<TATensor<std::complex<double>,4>> s;


    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma0_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma94_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma2_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma3_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma4_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma5_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma6_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma7_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma9_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma107_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma12_();
    std::shared_ptr<FutureTATensor<std::complex<double>,2>> Gamma14_();
    std::shared_ptr<FutureTATensor<std::complex<double>,2>> Gamma16_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma22_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma28_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma29_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma31_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma32_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma34_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma35_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma37_();
    std::shared_ptr<FutureTATensor<std::complex<double>,2>> Gamma38_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma51_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma56_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma57_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma58_();
    std::shared_ptr<FutureTATensor<std::complex<double>,6>> Gamma59_();
    std::shared_ptr<FutureTATensor_new<std::complex<double>,4>> Gamma60_();
    std::shared_ptr<FutureTATensor_new<std::complex<double>,0>> Gamma69_();
    std::shared_ptr<FutureTATensor<std::complex<double>,2>> Gamma81_();
    std::shared_ptr<FutureTATensor<std::complex<double>,4>> Gamma92_();
    std::shared_ptr<Queue> make_residualq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_sourceq(const bool reset = true, const bool diagonal = true);

  public:
    RelCASPT2(std::shared_ptr<const SMITH_Info<std::complex<double>>> ref);
    ~RelCASPT2() {}

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

