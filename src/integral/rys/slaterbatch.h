//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: slaterbatch.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_INTEGRAL_RYS_SLATERBATCH_H
#define __SRC_INTEGRAL_RYS_SLATERBATCH_H

#include <src/integral/rys/int2d.h>
#include <src/integral/rys/rysintegral.h>

namespace bagel {

#ifdef HAVE_LIBSLATER
class SlaterBatch : public RysInt {

  protected:
    double gamma_;
    bool yukawa_;

    /// buffer and intermediate storage
    double *bkup2_;
    std::vector<std::tuple<int, double, double>> indexpair23_;

    void perform_SVRR1();
    void perform_SVRR2();
    void perform_USVRR1();
    void perform_USVRR2();
    void perform_SVRR();
    void perform_USVRR();

    void root_weight(const int primsize_) override;
    void root1_direct();
    void root2_direct();

    void compute_ssss(const double) override;
    void allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) override;

  public:

    SlaterBatch(const std::array<std::shared_ptr<const Shell>,4>&, const double, const double gamma, const bool doyukawa = false,
                std::shared_ptr<StackMem> stack = nullptr);
    ~SlaterBatch();

    /// compute a batch of integrals
    void compute() override;

    constexpr static int nblocks() { return 1; }
};
#endif

}

#endif

