//
// BAGEL - Parallel electron correlation program.
// Filename: slaterbatch.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#ifndef __SRC_RYSINT_SLATERBATCH_H
#define __SRC_RYSINT_SLATERBATCH_H

#include <cassert>
#include <vector>
#include <src/rysint/int2d.h>
#include <src/rysint/rysint.h>
#include <memory>
#include <tuple>

namespace bagel {

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

  public:

    SlaterBatch(const std::array<std::shared_ptr<const Shell>,4>&, const double, const double gamma, const bool doyukawa = false);
    ~SlaterBatch();

    /// compute a batch of integrals
    void compute() override;

};

}

#endif

