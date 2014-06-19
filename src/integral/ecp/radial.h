//
// BAGEL - Parallel electron correlation program.
// Filename: radial.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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
    std::vector<double> integral() { return integral_; }
    double integral(const int ic = 0) const { return integral_.at(ic); }

    virtual std::vector<double> compute(const std::vector<double> r) = 0;

    void transform_Log(const int ngrid, const int m = 3);
    void transform_Ahlrichs(const int ngrid);
    void transform_Becke(const int ngrid);

    void GaussChebyshev2nd(const int ngrid);

};

}

#endif
