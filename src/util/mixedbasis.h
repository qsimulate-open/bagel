//
// BAGEL - Parallel electron correlation program.
// Filename: mixed_matrix1e.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_UTIL_MIXED_MATRIX1E_H
#define __SRC_UTIL_MIXED_MATRIX1E_H

#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <cassert>
#include <string>
#include <src/util/matrix.h>
#include <src/wfn/geometry.h>

namespace bagel {

template<typename T>
class MixedBasis : public Matrix {
  protected:

    void computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1) {
      // input = [b1, b0]
      const int dimb1 = input[0]->nbasis();
      const int dimb0 = input[1]->nbasis();
      T batch(input);
      batch.compute();

      copy_block(offsetb1, offsetb0, dimb1, dimb0, batch.data());
    }

  public:
    MixedBasis(const std::shared_ptr<const Geometry> g0, const std::shared_ptr<const Geometry> g1)
     : Matrix(g1->nbasis(), g0->nbasis()) {

      auto iter0 = g0->offsets().begin();
      for (auto& catom0 : g0->atoms()) {
        const std::vector<int> coffset0 = *iter0++;

        auto iter1 = g1->offsets().begin();
        for (auto& catom1 : g1->atoms()) {
          const std::vector<int> coffset1 = *iter1++;

          auto ioff0 = coffset0.begin();
          for (auto& b0 : catom0->shells()) {

            auto ioff1 = coffset1.begin();
            for (auto& b1 : catom1->shells()) {
              computebatch(std::array<std::shared_ptr<const Shell>,2>{{b1, b0}}, *ioff0, *ioff1++);
            }
            ++ioff0;
          }
        }
      }

    }
    ~MixedBasis() {};

    void print(const std::string in = "", const size_t size = 10) const {
      std::cout << "++++ " << in << " ++++" << std::endl;
      for (int i = 0; i != std::min(size,static_cast<size_t>(ndim_)); ++i) {
        for (int j = 0; j != std::min(size,static_cast<size_t>(mdim_)); ++j) {
          std::cout << std::fixed << std::setw(12) << std::setprecision(9) << data_[j*ndim_+i]  << " ";
        }
        std::cout << std::endl;
      }
    }

};

}

#endif
