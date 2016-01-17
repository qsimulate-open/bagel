//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: radial.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_INTEGRAL_ECP_RADIAL_H
#define __SRC_INTEGRAL_ECP_RADIAL_H

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <src/util/constants.h>
#include <src/util/timer.h>

namespace bagel {

class RadialInt {

  protected:
    int nc_;
    bool print_intermediate_;
    int max_iter_;
    double thresh_int_;
    std::vector<double> x_, w_, r_;
    std::vector<double> integral_;

  public:
    RadialInt(const int nc = 1, const bool print = false, const int max_iter = 100, const double thresh_int = PRIM_SCREEN_THRESH) :
      nc_(nc),
      print_intermediate_(print),
      max_iter_(max_iter),
      thresh_int_(thresh_int)
    {}

    ~RadialInt() {}

    void integrate();
    std::vector<double> integral() const { return integral_; }
    double integral(const int ic = 0) const { return integral_.at(ic); }

    virtual std::vector<double> compute(const std::vector<double> r) = 0;

    void transform_Log(const int ngrid, const int m = 3);
    void transform_Ahlrichs(const int ngrid);
    void transform_Becke(const int ngrid);

    void GaussChebyshev2nd(const int ngrid);

};

}

#endif
