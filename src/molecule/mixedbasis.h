//
// BAGEL - Parallel electron correlation program.
// Filename: mixedbasis.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_MOLECULE_MIXEDBASIS_H
#define __SRC_MOLECULE_MIXEDBASIS_H

#include <src/molecule/molecule.h>

namespace bagel {

// variadic template
template<typename T, typename... Value>
class MixedBasis : public Matrix {
  protected:

    void computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, Value... tail) {
      // input = [b1, b0]
      const int dimb1 = input[0]->nbasis();
      const int dimb0 = input[1]->nbasis();
      T batch(input, tail...);
      batch.compute();

      copy_block(offsetb1, offsetb0, dimb1, dimb0, batch.data());
    }

  public:
    MixedBasis(const std::shared_ptr<const Molecule> g0, const std::shared_ptr<const Geometry> g1, Value... tail)
     : Matrix(g1->nbasis(), g0->nbasis()) {
      size_t off0 = 0;
      for (auto& catom0 : g0->atoms()) {
        for (auto& b0 : catom0->shells()) {
          size_t off1 = 0;
          for (auto& catom1 : g1->atoms())
            for (auto& b1 : catom1->shells()) {
              computebatch(std::array<std::shared_ptr<const Shell>,2>{{b1, b0}}, off0, off1, tail...);
              off1 += b1->nbasis();
            }
          off0 += b0->nbasis();
        }
      }
    }

    void print(const std::string in = "", const int size = 10) const {
      std::cout << "++++ " << in << " ++++" << std::endl;
      for (int i = 0; i != std::min(size, ndim()); ++i) {
        for (int j = 0; j != std::min(size, mdim()); ++j) {
          std::cout << std::fixed << std::setw(12) << std::setprecision(9) << element(i, j) << " ";
        }
        std::cout << std::endl;
      }
    }

};

}

#endif
