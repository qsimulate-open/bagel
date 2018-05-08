//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI.h
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
#include <src/grad/nacmtype.h>

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

    std::vector<std::shared_ptr<MultiTensor>> t2all_;
    std::vector<std::shared_ptr<MultiTensor>> rall_;
    std::vector<std::shared_ptr<MultiTensor>> sall_;
    std::vector<std::shared_ptr<MultiTensor>> nall_;
    void diagonal(std::shared_ptr<Tensor> r, std::shared_ptr<const Tensor> t) const;


    std::shared_ptr<FutureTensor> Gamma0_();
    std::shared_ptr<FutureTensor> Gamma1_();
    std::shared_ptr<FutureTensor> Gamma2_();
    std::shared_ptr<FutureTensor> Gamma58_();
    std::shared_ptr<FutureTensor> Gamma59_();
    std::shared_ptr<FutureTensor> Gamma60_();
    std::shared_ptr<FutureTensor> Gamma63_();
    std::shared_ptr<FutureTensor> Gamma64_();
    std::shared_ptr<FutureTensor> Gamma65_();
    std::shared_ptr<FutureTensor> Gamma66_();
    std::shared_ptr<FutureTensor> Gamma67_();
    std::shared_ptr<FutureTensor> Gamma3_();
    std::shared_ptr<FutureTensor> Gamma4_();
    std::shared_ptr<FutureTensor> Gamma5_();
    std::shared_ptr<FutureTensor> Gamma74_();
    std::shared_ptr<FutureTensor> Gamma75_();
    std::shared_ptr<FutureTensor> Gamma77_();
    std::shared_ptr<FutureTensor> Gamma78_();
    std::shared_ptr<FutureTensor> Gamma79_();
    std::shared_ptr<FutureTensor> Gamma81_();
    std::shared_ptr<FutureTensor> Gamma84_();
    std::shared_ptr<FutureTensor> Gamma86_();
    std::shared_ptr<FutureTensor> Gamma92_();
    std::shared_ptr<FutureTensor> Gamma95_();
    std::shared_ptr<FutureTensor> Gamma98_();
    std::shared_ptr<FutureTensor> Gamma414_();
    std::shared_ptr<FutureTensor> Gamma415_();
    std::shared_ptr<FutureTensor> Gamma9_();
    std::shared_ptr<FutureTensor> Gamma11_();
    std::shared_ptr<FutureTensor> Gamma160_();
    std::shared_ptr<FutureTensor> Gamma99_();
    std::shared_ptr<FutureTensor> Gamma105_();
    std::shared_ptr<FutureTensor> Gamma110_();
    std::shared_ptr<FutureTensor> Gamma128_();
    std::shared_ptr<FutureTensor> Gamma151_();
    std::shared_ptr<FutureTensor> Gamma159_();
    std::shared_ptr<FutureTensor> Gamma174_();
    std::shared_ptr<FutureTensor> Gamma416_();
    std::shared_ptr<FutureTensor> Gamma418_();
    std::shared_ptr<FutureTensor> Gamma24_();
    std::shared_ptr<FutureTensor> Gamma26_();
    std::shared_ptr<FutureTensor> Gamma27_();
    std::shared_ptr<FutureTensor> Gamma179_();
    std::shared_ptr<FutureTensor> Gamma183_();
    std::shared_ptr<FutureTensor> Gamma193_();
    std::shared_ptr<FutureTensor> Gamma203_();
    std::shared_ptr<FutureTensor> Gamma204_();
    std::shared_ptr<FutureTensor> Gamma229_();
    std::shared_ptr<FutureTensor> Gamma420_();
    std::shared_ptr<FutureTensor> Gamma421_();
    std::shared_ptr<FutureTensor> Gamma31_();
    std::shared_ptr<FutureTensor> Gamma32_();
    std::shared_ptr<FutureTensor> Gamma33_();
    std::shared_ptr<FutureTensor> Gamma230_();
    std::shared_ptr<FutureTensor> Gamma231_();
    std::shared_ptr<FutureTensor> Gamma232_();
    std::shared_ptr<FutureTensor> Gamma233_();
    std::shared_ptr<FutureTensor> Gamma236_();
    std::shared_ptr<FutureTensor> Gamma237_();
    std::shared_ptr<FutureTensor> Gamma238_();
    std::shared_ptr<FutureTensor> Gamma245_();
    std::shared_ptr<FutureTensor> Gamma246_();
    std::shared_ptr<FutureTensor> Gamma253_();
    std::shared_ptr<FutureTensor> Gamma422_();
    std::shared_ptr<FutureTensor> Gamma423_();
    std::shared_ptr<FutureTensor> Gamma317_();
    std::shared_ptr<FutureTensor> Gamma318_();
    std::shared_ptr<FutureTensor> Gamma335_();
    std::shared_ptr<FutureTensor> Gamma336_();
    std::shared_ptr<FutureTensor> Gamma368_();
    std::shared_ptr<FutureTensor> Gamma363_();
    std::shared_ptr<FutureTensor> Gamma391_();
    std::shared_ptr<FutureTensor> Gamma428_();
    std::shared_ptr<FutureTensor> Gamma430_();
    std::shared_ptr<FutureTensor> Gamma396_();
    std::shared_ptr<FutureTensor> Gamma397_();
    std::shared_ptr<FutureTensor> Gamma409_();
    std::shared_ptr<FutureTensor> Gamma412_();
    std::shared_ptr<FutureTensor> Gamma413_();
    std::shared_ptr<FutureTensor> Gamma424_();
    std::shared_ptr<FutureTensor> Gamma426_();
    std::shared_ptr<FutureTensor> Gamma432_();
    std::shared_ptr<FutureTensor> Gamma433_();

    std::shared_ptr<Queue> make_residualq(const bool reset = true, const bool diagonal = true);
    void make_residualq1(std::shared_ptr<Queue> out, std::shared_ptr<Task> task, const bool diagonal);
    void make_residualq2(std::shared_ptr<Queue> out, std::shared_ptr<Task> task, const bool diagonal);
    void make_residualq3(std::shared_ptr<Queue> out, std::shared_ptr<Task> task, const bool diagonal);
    void make_residualq4(std::shared_ptr<Queue> out, std::shared_ptr<Task> task, const bool diagonal);
    void make_residualq5(std::shared_ptr<Queue> out, std::shared_ptr<Task> task, const bool diagonal);
    void make_residualq6(std::shared_ptr<Queue> out, std::shared_ptr<Task> task, const bool diagonal);

    std::shared_ptr<Queue> make_sourceq(const bool reset = true, const bool diagonal = true);
    std::shared_ptr<Queue> make_normq(const bool reset = true, const bool diagonal = true);

  public:
    RelMRCI(std::shared_ptr<const SMITH_Info<std::complex<double>>> ref);
    ~RelMRCI() {}

    void solve();
    void solve_gradient(const int targetJ, const int targetI, std::shared_ptr<const NacmType> nacmtype = std::make_shared<const NacmType>(), const bool nocider = false);

    double accumulate(std::shared_ptr<Queue> queue) {
      double sum = 0.0;
      while (!queue->done())
        sum += queue->next_compute()->target();
      mpi__->allreduce(&sum, 1);
      return sum;
    }

};

}
}
}
#endif

