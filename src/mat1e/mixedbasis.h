//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mixedbasis.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_MAT1E_MIXEDBASIS_H
#define __SRC_MAT1E_MIXEDBASIS_H

#include <src/molecule/molecule.h>

namespace bagel {

// variadic template
template<typename TBatch, typename MatType=Matrix, typename... Value>
class MixedBasis : public MatType {
  protected:

    void computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, Value... tail) {
      // input = [b1, b0]
      const int dimb1 = input[0]->nbasis();
      const int dimb0 = input[1]->nbasis();
      TBatch batch(input, tail...);
      batch.compute();

      this->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch.data());
    }

  public:
    MixedBasis(std::shared_ptr<const Molecule> g0, std::shared_ptr<const Molecule> g1, Value... tail)
     : MatType(g1->nbasis(), g0->nbasis()) {
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

    void print(const std::string in = "", const size_t size = 10) const {
      std::cout << "++++ " << in << " ++++" << std::endl;
      for (int i = 0; i != std::min(size, this->ndim()); ++i) {
        for (int j = 0; j != std::min(size, this->mdim()); ++j) {
          std::cout << std::fixed << std::setw(12) << std::setprecision(9) << this->element(i, j) << " ";
        }
        std::cout << std::endl;
      }
    }

};


template<typename TBatch, typename MatType=Matrix, typename... Value>
class MixedBasisArray {
  protected:

    std::array<std::shared_ptr<MatType>,TBatch::Nblocks()> matrices_;

    void computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, Value... tail) {
      // input = [b1, b0]
      const int dimb1 = input[0]->nbasis();
      const int dimb0 = input[1]->nbasis();
      TBatch batch(input, tail...);
      batch.compute();

      for (int i = 0; i != TBatch::Nblocks(); ++i)
        matrices_[i]->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch[i]);
    }

  public:
    MixedBasisArray(std::shared_ptr<const Molecule> g0, std::shared_ptr<const Molecule> g1, Value... tail) {
      for(int i = 0; i < TBatch::Nblocks(); ++i)
        matrices_[i] = std::make_shared<MatType>(g1->nbasis(), g0->nbasis());

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

    constexpr static int Nblocks() { return TBatch::Nblocks(); }

    std::shared_ptr<const MatType> data(const int i) const { return matrices_[i]; }

    void print(const std::string in = "", const int size = 10) const {
      for (int k = 0; k != TBatch::Nblocks(); ++k) {
        std::cout << "++++ " << in << " ++++ Block " << k << " ++++" << std::endl;
        for (int i = 0; i != std::min(size, matrices_[k]->ndim()); ++i) {
          for (int j = 0; j != std::min(size, matrices_[k]->mdim()); ++j) {
            std::cout << std::fixed << std::setw(12) << std::setprecision(9) << matrices_[k]->element(i, j) << " ";
          }
          std::cout << std::endl;
        }
      }
    }

};

}

#endif
