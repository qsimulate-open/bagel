//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: smallints1e_general.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

// A specialized small component integral class for when we need nonstandard combinations (such as z^z - x^x - y^y and x^z + z^x)
// Presently only used for spin expectation values
// Very similar to standard SmallInts1e - TODO Eliminate redundancy


#ifndef __SRC_INTEGRAL_SMALLINTS1E_GENERAL_H
#define __SRC_INTEGRAL_SMALLINTS1E_GENERAL_H

#include <src/util/math/matrix.h>

// computes (sigma p)Vnuc(sigma p), and returns 9 blocks of data

namespace bagel {

template<typename Batch>
class SmallInts1e_General {
  protected:
    std::array<std::shared_ptr<Matrix>,9*Batch::Nblocks()> data_;

    const std::array<std::shared_ptr<const Shell>,2> shells_;

    const size_t size_block_;

    template<typename Value>
    std::array<std::shared_ptr<Matrix>,Batch::Nblocks()> int_compute(const Value&) const {
    }

    void transform(const std::array<std::shared_ptr<Matrix>,Batch::Nblocks()>& unc) {
      constexpr int N = Batch::Nblocks();

      for (int n = 0; n != N; ++n) {
        std::array<std::shared_ptr<Matrix>,3> ints;
        for (int i = 0; i != 3; ++i)
          ints[i] = std::make_shared<Matrix>(*shells_[0]->small(i) % *unc[n]);

        std::array<int,3> f = {{2,3,1}};
        std::array<int,3> b = {{3,1,2}};

        // 0) x^x
        // 1) y^y
        // 2) z^z

        // 3) x^y
        // 4) y^z
        // 5) z^x

        // 6) y^x
        // 7) z^y
        // 8) x^z

        for (int i=0; i!=3; ++i) {
          *data_[9*n+i]      = *ints[i]      * *shells_[1]->small(i);
          *data_[9*n+b[i]+2] = *ints[b[i]-1] * *shells_[1]->small(i);
          *data_[9*n+i+6]    = *ints[f[i]-1] * *shells_[1]->small(i);
        }
      }
    }

  public:
    SmallInts1e_General(std::array<std::shared_ptr<const Shell>,2> info)
      : shells_(info), size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()) {

      for (int i = 0; i != Nblocks(); ++i)
        data_[i] = std::make_shared<Matrix>(shells_[0]->nbasis(), shells_[1]->nbasis(), true);
    }

    void compute() { compute<void*>(nullptr); }

    template<typename Value>
    void compute(const Value&) {
      static_assert(std::is_same<Value, void*>::value, "SmallInts1e_General::compute called illegally");
      const int a0size_inc = shells_[0]->nbasis_aux_increment();
      const int a1size_inc = shells_[1]->nbasis_aux_increment();
      const int a0size_dec = shells_[0]->nbasis_aux_decrement();
      const int a1size_dec = shells_[1]->nbasis_aux_decrement();
      const int a0 = a0size_inc + a0size_dec;
      const int a1 = a1size_inc + a1size_dec;

      constexpr int N = Batch::Nblocks();

      std::array<std::shared_ptr<Matrix>,N> unc;
      for (int n = 0; n != N; ++n)
        unc[n] = std::make_shared<Matrix>(a0, a1, true);

      {
        auto uncc = std::make_shared<Batch>(std::array<std::shared_ptr<const Shell>,2>{{shells_[0]->aux_increment(), shells_[1]->aux_increment()}});
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(0, 0, a0size_inc, a1size_inc, uncc->data(n));
      }
      if (shells_[0]->aux_decrement() && shells_[1]->aux_decrement()) {
        auto uncc = std::make_shared<Batch>(std::array<std::shared_ptr<const Shell>,2>{{shells_[0]->aux_decrement(), shells_[1]->aux_decrement()}});
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(a0size_inc, a1size_inc, a0size_dec, a1size_dec, uncc->data(n));
      }
      if (shells_[0]->aux_decrement()) {
        auto uncc = std::make_shared<Batch>(std::array<std::shared_ptr<const Shell>,2>{{shells_[0]->aux_decrement(), shells_[1]->aux_increment()}});
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(a0size_inc, 0, a0size_dec, a1size_inc, uncc->data(n));
      }
      if (shells_[1]->aux_decrement()) {
        auto uncc = std::make_shared<Batch>(std::array<std::shared_ptr<const Shell>,2>{{shells_[0]->aux_increment(), shells_[1]->aux_decrement()}});
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(0, a1size_inc, a0size_inc, a1size_dec, uncc->data(n));
      }

      transform(unc);
    }

    std::shared_ptr<Matrix> operator[](const int i) { return data_[i]; }

    size_t size_block() const { return size_block_; }
    constexpr static int Nblocks() { return Batch::Nblocks() * 9; }

};

}

#endif
