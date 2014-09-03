//
// BAGEL - Parallel electron correlation program.
// Filename: soangularbatch.h
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


#ifndef __SRC_INTEGRAL_ECP_SOBATCH_H
#define __SRC_INTEGRAL_ECP_SOBATCH_H

#include <src/molecule/atom.h>
#include <src/integral/ecp/sphharmonics.h>
#include <src/integral/ecp/radial.h>
#include <map>

namespace bagel {

class SOBatch : public RadialInt {
  protected:

    std::array<std::shared_ptr<const Shell>,2> basisinfo_;
    std::shared_ptr<const SOECP> so_;
    int cont0_, cont1_;
    std::array<int, 3> ang0_, ang1_;
    std::array<double, 3> AB_, CB_;
    double dAB_, dCB_;
    std::vector<double> c0_, c1_;
    std::vector<std::map<int, std::array<int, 3>>> map_;

    int l0_, l1_;
    std::vector<std::vector<double>> zAB_, zCB_;
    std::vector<std::vector<std::tuple<int, int, int, double>>> fm0lm1_; // li, m0, m1, fmm

    void map_angular_number();
    std::complex<double> theta(const int m) const;

    std::array<double, 3> fm0lm1(const int l, const int m0, const int m1) const;
    std::vector<double> project(const int l, const std::vector<double> r);
    double angularA(const int h, const int ld, const std::vector<double> usp);
    double angularC(const int h, const int ld, const std::vector<double> usp);

  public:
    SOBatch(const std::shared_ptr<const SOECP> _so, const std::array<std::shared_ptr<const Shell>,2>& _info,
                 const double contA, const double contC, const std::array<int, 3> angA, const std::array<int, 3> angC,
                 const bool print = false, const int max_iter = 100, const double thresh_int = PRIM_SCREEN_THRESH);

    ~SOBatch() {}

    std::vector<double> compute(const std::vector<double> r) override;

    void init();
    void print() const;

};

}

#endif
