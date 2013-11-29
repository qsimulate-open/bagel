
// BAGEL - Parallel electron correlation program.
// Filename: moldenin.h
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __SRC_IO_MOLDENIN_H
#define __SRC_IO_MOLDENIN_H

#include <src/io/moldenio.h>

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
        MoldenIn& operator>> (std::shared_ptr<Coeff>& coeff);
  };
}

#endif
