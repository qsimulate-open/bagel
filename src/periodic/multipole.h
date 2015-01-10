//
// BAGEL - Parallel electron correlation program.
// Filename: multipole.h
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


#ifndef __SRC_PERIODIC_MULTIPOLE_H
#define __SRC_PERIODIC_MULTIPOLE_H

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <src/math/legendre.h>
#include <src/util/constants.h>

namespace bagel {

class Multipole {
  protected:
    std::array<double, 3> centre_;
    int lmax_;

    int num_multipoles_;
    std::vector<std::complex<double>> multipole_;
    void compute_multipoles();

  public:
    Multipole() { }
    Multipole(const std::array<double, 3> c, const int lmax = ANG_HRR_END);
    ~Multipole() { }

    std::array<double, 3> centre() const { return centre_; }
    double centre(const int i) const { return centre_[i]; }

    int num_multipoles() const { return num_multipoles_; }

    std::vector<std::complex<double>> multipoles() { return multipole_; }
    std::complex<double> multipole(const int i) const { assert(i < num_multipoles_); return multipole_[i]; }
    std::complex<double> multipole(const int l, const int m) const {
      assert (l <= lmax_ && abs(m) <= l); return multipole_[l * l + l + m];
    }

    std::vector<std::complex<double>> multipoles(const int l) {
      assert (l <= lmax_);
      std::vector<std::complex<double>> out(2 * l + 1);
      const int i0 = (l + 1) * (l + 1);
      for (int i = 0; i != 2 * l + 1; ++i) out[i] = multipole_[i + i0];

      return out;
    }

//    std::shared_ptr<const Multipole> translate(std::array<double, 3> displ);

    void print_multipoles() const;
};

}

#endif
