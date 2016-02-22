//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: smallints1e_london.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __SRC_INTEGRAL_SMALLINTS1E_LONDON_H
#define __SRC_INTEGRAL_SMALLINTS1E_LONDON_H

#include <src/integral/comprys/complexeribatch.h>
#include <src/molecule/molecule.h>
#include <src/util/unpack.h>

// computes (sigma pi)Operator(sigma pi), and returns 4 blocks of data

namespace bagel {

template<typename Batch, typename ...Args>
class SmallInts1e_London {
  protected:
    std::array<std::shared_ptr<ZMatrix>,4*Batch::Nblocks()> data_;

    const std::tuple<Args...> args_;
    const std::array<std::shared_ptr<const Shell>,2> shells_;

    const size_t size_block_;

    template<typename Value>
    std::array<std::shared_ptr<ZMatrix>,Batch::Nblocks()> int_compute(const Value&) const {
    }

    void transform(const std::array<std::shared_ptr<ZMatrix>,Batch::Nblocks()>& unc) {
      constexpr int N = Batch::Nblocks();

      for (int n = 0; n != N; ++n) {
        std::array<std::shared_ptr<ZMatrix>,3> ints;
        for (int i = 0; i != 3; ++i) {
          ints[i] = std::make_shared<ZMatrix>(*shells_[0]->zsmall(i) % *unc[n]);
        }

        std::array<int,3> f = {{2,3,1}};
        std::array<int,3> b = {{3,1,2}};

        // 0) x^x + y^y + z^z
        // 1) x^y - y^x
        // 2) y^z - z^y
        // 3) z^x - x^z

        // -1 because <m|p|n>^dagger = -<n|p|m>  (can be proven by integration by part)
        for (int i = 0; i != 3; ++i) {
          *data_[4*n+0]    += *ints[i]      * *shells_[1]->zsmall(i);
          *data_[4*n+b[i]] += *ints[b[i]-1] * *shells_[1]->zsmall(i);
          *data_[4*n+i+1]  -= *ints[f[i]-1] * *shells_[1]->zsmall(i);
        }
      }
    }

  public:
    SmallInts1e_London(std::array<std::shared_ptr<const Shell>,2> info, Args... args)
      : args_(args...), shells_(info), size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()) {

      for (int i = 0; i != Nblocks(); ++i)
        data_[i] = std::make_shared<ZMatrix>(shells_[0]->nbasis(), shells_[1]->nbasis(), true);
    }

    void compute() { compute<void*>(nullptr); }

    // Unpack variadic template arguments here
    template <int ...S>
    std::shared_ptr<Batch> get_batch(std::shared_ptr<const Shell> a0, std::shared_ptr<const Shell> a1, seq<S...>) {
      auto out = std::make_shared<Batch>(std::array<std::shared_ptr<const Shell>,2>{{a0, a1}}, std::get<S>(args_) ...);
      return out;
    }

    template<typename Value>
    void compute(const Value&) {
      static_assert(std::is_same<Value, void*>::value, "SmallInts1e_London::compute called illegally");
      const int a0size_inc = shells_[0]->nbasis_aux_increment();
      const int a1size_inc = shells_[1]->nbasis_aux_increment();
      const int a0size_dec = shells_[0]->nbasis_aux_decrement();
      const int a1size_dec = shells_[1]->nbasis_aux_decrement();
      const int a0size_same = shells_[0]->nbasis_aux_same();
      const int a1size_same = shells_[1]->nbasis_aux_same();
      const int a0 = a0size_inc + a0size_dec + a0size_same;
      const int a1 = a1size_inc + a1size_dec + a1size_same;

      constexpr int N = Batch::Nblocks();

      std::array<std::shared_ptr<ZMatrix>,N> unc;
      for (int n = 0; n != N; ++n)
        unc[n] = std::make_shared<ZMatrix>(a0, a1, true);

      {
        auto uncc = get_batch(shells_[0]->aux_increment(), shells_[1]->aux_increment(), typename gens<sizeof...(Args)>::type());
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(0, 0, a0size_inc, a1size_inc, uncc->data(n));
      }
      if (shells_[0]->aux_decrement() && shells_[1]->aux_decrement()) {
        auto uncc = get_batch(shells_[0]->aux_decrement(), shells_[1]->aux_decrement(), typename gens<sizeof...(Args)>::type());
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(a0size_inc, a1size_inc, a0size_dec, a1size_dec, uncc->data(n));
      }
      if (shells_[0]->aux_decrement()) {
        auto uncc = get_batch(shells_[0]->aux_decrement(), shells_[1]->aux_increment(), typename gens<sizeof...(Args)>::type());
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(a0size_inc, 0, a0size_dec, a1size_inc, uncc->data(n));
      }
      if (shells_[1]->aux_decrement()) {
        auto uncc = get_batch(shells_[0]->aux_increment(), shells_[1]->aux_decrement(), typename gens<sizeof...(Args)>::type());
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(0, a1size_inc, a0size_inc, a1size_dec, uncc->data(n));
      }

      // Unchanged angular momentum (common origin only)
      if (shells_[0]->aux_same()) {
        assert(shells_[1]->aux_same());
        const size_t a0size_id = a0size_inc + a0size_dec;
        const size_t a1size_id = a1size_inc + a1size_dec;
        {
          auto uncc = get_batch(shells_[0]->aux_increment(), shells_[1]->aux_same(), typename gens<sizeof...(Args)>::type());
          uncc->compute();
          for (int n = 0; n != N; ++n)
            unc[n]->copy_block(0, a1size_id, a0size_inc, a1size_same, uncc->data(n));
        }
        {
          auto uncc = get_batch(shells_[0]->aux_same(), shells_[1]->aux_increment(), typename gens<sizeof...(Args)>::type());
          uncc->compute();
          for (int n = 0; n != N; ++n)
            unc[n]->copy_block(a0size_id, 0, a0size_same, a1size_inc, uncc->data(n));
        }
        {
          auto uncc = get_batch(shells_[0]->aux_same(), shells_[1]->aux_same(), typename gens<sizeof...(Args)>::type());
          uncc->compute();
          for (int n = 0; n != N; ++n)
            unc[n]->copy_block(a0size_id, a1size_id, a0size_same, a1size_same, uncc->data(n));
        }
        if (shells_[0]->aux_decrement()) {
          auto uncc = get_batch(shells_[0]->aux_decrement(), shells_[1]->aux_same(), typename gens<sizeof...(Args)>::type());
          uncc->compute();
          for (int n = 0; n != N; ++n)
            unc[n]->copy_block(a0size_inc, a1size_id, a0size_dec, a1size_same, uncc->data(n));
        }
        if (shells_[1]->aux_decrement()) {
          auto uncc = get_batch(shells_[0]->aux_same(), shells_[1]->aux_decrement(), typename gens<sizeof...(Args)>::type());
          uncc->compute();
          for (int n = 0; n != N; ++n)
            unc[n]->copy_block(a0size_id, a1size_inc, a0size_same, a1size_dec, uncc->data(n));
        }
      } else assert(!shells_[1]->aux_same());

      transform(unc);
    }

    std::shared_ptr<ZMatrix> operator[](const int i) { return data_[i]; }

    size_t size_block() const { return size_block_; }
    constexpr static int Nblocks() { return Batch::Nblocks() * 4; }

};


// in the case of finite nucleus
template<>
template<typename Value>
void SmallInts1e_London<ComplexERIBatch, std::shared_ptr<const Molecule>>::compute(const Value& nshells) {
  const int a0size_inc = shells_[0]->nbasis_aux_increment();
  const int a1size_inc = shells_[1]->nbasis_aux_increment();
  const int a0size_dec = shells_[0]->nbasis_aux_decrement();
  const int a1size_dec = shells_[1]->nbasis_aux_decrement();
  const int a0size_same = shells_[0]->nbasis_aux_same();
  const int a1size_same = shells_[1]->nbasis_aux_same();
  const int a0 = a0size_inc + a0size_dec + a0size_same;
  const int a1 = a1size_inc + a1size_dec + a1size_same;

  auto dummy = std::make_shared<const Shell>(shells_[0]->spherical());

  constexpr int N = ComplexERIBatch::Nblocks();
  static_assert(N == 1, "ComplexERIBatch::Nblocks() should be 1");

  std::array<std::shared_ptr<ZMatrix>,N> unc;
  for (int n = 0; n != N; ++n)
    unc[n] = std::make_shared<ZMatrix>(a0, a1, true);

  for (auto& nshell : nshells) {
    {
      auto uncc = std::make_shared<ComplexERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_increment(), shells_[1]->aux_increment()}}, 2.0);
      uncc->compute();
      for (int n = 0; n != N; ++n)
        unc[n]->add_block(1.0, 0, 0, a0size_inc, a1size_inc, uncc->data(n));
    }
    if (shells_[0]->aux_decrement() && shells_[1]->aux_decrement()) {
      auto uncc = std::make_shared<ComplexERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_decrement(), shells_[1]->aux_decrement()}}, 2.0);
      uncc->compute();
      for (int n = 0; n != N; ++n)
        unc[n]->add_block(1.0, a0size_inc, a1size_inc, a0size_dec, a1size_dec, uncc->data(n));
    }
    if (shells_[0]->aux_decrement()) {
      auto uncc = std::make_shared<ComplexERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_decrement(), shells_[1]->aux_increment()}}, 2.0);
      uncc->compute();
      for (int n = 0; n != N; ++n)
        unc[n]->add_block(1.0, a0size_inc, 0, a0size_dec, a1size_inc, uncc->data(n));
    }
    if (shells_[1]->aux_decrement()) {
      auto uncc = std::make_shared<ComplexERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_increment(), shells_[1]->aux_decrement()}}, 2.0);
      uncc->compute();
      for (int n = 0; n != N; ++n)
        unc[n]->add_block(1.0, 0, a1size_inc, a0size_inc, a1size_dec, uncc->data(n));
    }

    // Unchanged angular momentum (common origin only)
    if (shells_[0]->aux_same()) {
      assert(shells_[1]->aux_same());
      const size_t a0size_id = a0size_inc + a0size_dec;
      const size_t a1size_id = a1size_inc + a1size_dec;
      {
        auto uncc = std::make_shared<ComplexERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_increment(), shells_[1]->aux_same()}}, 2.0);
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->add_block(1.0, 0, a1size_id, a0size_inc, a1size_same, uncc->data(n));
      }
      {
        auto uncc = std::make_shared<ComplexERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_same(), shells_[1]->aux_increment()}}, 2.0);
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->add_block(1.0, a0size_id, 0, a0size_same, a1size_inc, uncc->data(n));
      }
      {
        auto uncc = std::make_shared<ComplexERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_same(), shells_[1]->aux_same()}}, 2.0);
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->add_block(1.0, a0size_id, a1size_id, a0size_same, a1size_same, uncc->data(n));
      }
      if (shells_[0]->aux_decrement()) {
        auto uncc = std::make_shared<ComplexERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_decrement(), shells_[1]->aux_same()}}, 2.0);
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->add_block(1.0, a0size_inc, a1size_id, a0size_dec, a1size_same, uncc->data(n));
      }
      if (shells_[1]->aux_decrement()) {
        auto uncc = std::make_shared<ComplexERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_same(), shells_[1]->aux_decrement()}}, 2.0);
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->add_block(1.0, a0size_id, a1size_inc, a0size_same, a1size_dec, uncc->data(n));
       }
    } else assert(!shells_[1]->aux_same());
  }

  transform(unc);
}

}

#endif
