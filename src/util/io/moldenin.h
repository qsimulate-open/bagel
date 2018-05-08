//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moldenin.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __SRC_IO_MOLDENIN_H
#define __SRC_IO_MOLDENIN_H

#include <src/util/io/moldenio.h>

namespace bagel {

class MoldenIn : public MoldenIO {
  protected:
    bool is_spherical_;
    bool cartesian_;

    std::vector<std::shared_ptr<const Atom>> atoms_;

    std::vector<std::vector<double>> mo_coefficients_;

    std::vector<std::vector<std::vector<std::pair<int, double>>>> lmtuv_;
    std::vector<int> gto_order_;
    std::vector<std::vector<int>> shell_orders_;

    void compute_transforms();
    std::vector<double> transform_cart(std::vector<double> in, int ang_l);

  public:
    MoldenIn(const std::string filename, const bool is_spherical = true);

    void read();

    MoldenIn& operator>> (std::vector<std::shared_ptr<const Atom>>& atoms_);
    MoldenIn& operator>> (std::tuple<std::shared_ptr<Coeff>, std::shared_ptr<const Geometry>>);
    bool has_mo() const { return !mo_coefficients_.empty(); }
};

}

#endif
