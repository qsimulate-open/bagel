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

class CASPT2 : public SpinFreeMethod<double> {
  protected:
    std::shared_ptr<TATensor<double,4>> t2;
    std::shared_ptr<TATensor<double,4>> r;
    std::shared_ptr<TATensor<double,2>> den1;
    std::shared_ptr<TATensor<double,2>> den2;
    std::shared_ptr<TATensor<double,4>> Den1;
    double correlated_norm_;
    std::shared_ptr<TATensor<double,1>> deci;

    void diagonal(std::shared_ptr<TATensor<double,4>> r, std::shared_ptr<const TATensor<double,4>> t) const;

    std::shared_ptr<FutureTATensor<double,4>> Gamma0_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma92_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma2_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma3_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma4_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma5_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma6_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma7_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma9_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma105_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma12_();
    std::shared_ptr<FutureTATensor<double,2>> Gamma14_();
    std::shared_ptr<FutureTATensor<double,2>> Gamma16_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma22_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma28_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma29_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma31_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma32_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma34_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma35_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma37_();
    std::shared_ptr<FutureTATensor<double,2>> Gamma38_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma51_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma56_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma57_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma58_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma59_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma60_();
    std::shared_ptr<FutureTATensor<double,2>> Gamma79_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma90_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma160_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma191_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma194_();
    std::shared_ptr<FutureTATensor<double,6>> Gamma252_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma165_();
    std::shared_ptr<FutureTATensor<double,8>> Gamma218_();
    std::shared_ptr<FutureTATensor<double,4>> Gamma174_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma270_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma271_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma272_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma273_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma274_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma275_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma276_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma277_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma279_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma282_();
    std::shared_ptr<FutureTATensor<double,3>> Gamma284_();
    std::shared_ptr<FutureTATensor<double,3>> Gamma286_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma292_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma298_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma299_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma301_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma302_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma304_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma305_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma307_();
    std::shared_ptr<FutureTATensor<double,3>> Gamma308_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma321_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma326_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma327_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma328_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma329_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma330_();
    std::shared_ptr<FutureTATensor<double,1>> Gamma339_();
    std::shared_ptr<FutureTATensor<double,3>> Gamma351_();
    std::shared_ptr<FutureTATensor<double,5>> Gamma362_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma377_();
    std::shared_ptr<FutureTATensor<double,7>> Gamma395_();
    std::shared_ptr<Queue> make_residualq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_energyq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_corrq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_densityq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_density1q(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_density2q(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_deciq(const bool reset = true, const bool diagonal = true);

  public:
    CASPT2(std::shared_ptr<const SMITH_Info<double>> ref);
    ~CASPT2() {}

    void solve();
    void solve_deriv();

    double accumulate(std::shared_ptr<Queue> queue) {
      double sum = 0.0;
      while (!queue->done())
        sum += queue->next_compute()->target();
      return sum;
    }

    std::shared_ptr<const Matrix> rdm11() const { const Tensor d(*den1); return d.matrix(); }
    std::shared_ptr<const Matrix> rdm12() const { const Tensor d(*den2); return d.matrix(); }
    std::shared_ptr<const Matrix> rdm21() const { const Tensor d(*Den1); return d.matrix2(); }

    double correlated_norm() const { return correlated_norm_; }

    std::shared_ptr<const Civec> ci_deriv(std::shared_ptr<const Determinants> det) const { const Tensor d(*deci); return d.civec(det); }

};

}
}
}
#endif

