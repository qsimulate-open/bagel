//
// BAGEL - Parallel electron correlation program.
// Filename: smallints1e.h
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


#ifndef __SRC_INTEGRAL_SMALLINTS1E_H
#define __SRC_INTEGRAL_SMALLINTS1E_H

#include <src/molecule/molecule.h>
#include <src/integral/rys/eribatch.h>
#include <src/integral/libint/libint.h>

// computes (sigma p)Vnuc(sigma p), and returns 4 blocks of data

namespace bagel {

template<typename Batch>
class SmallInts1e {
  protected:
    std::array<std::shared_ptr<Matrix>,4*Batch::Nblocks()> data_;

    const std::shared_ptr<const Molecule> mol_;
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

        // 0) x^x + y^y + z^z
        // 1) x^y - y^x
        // 2) y^z - z^y
        // 3) z^x - x^z

        // -1 because <m|p|n>^dagger = -<n|p|m>  (can be proven by integration by part)
        for (int i = 0; i != 3; ++i) {
          *data_[4*n+0]    += *ints[i]      * *shells_[1]->small(i);
          *data_[4*n+b[i]] += *ints[b[i]-1] * *shells_[1]->small(i);
          *data_[4*n+i+1]  -= *ints[f[i]-1] * *shells_[1]->small(i);
        }
      }
    }

  public:
    SmallInts1e(std::array<std::shared_ptr<const Shell>,2> info, std::shared_ptr<const Molecule> mol)
      : mol_(mol), shells_(info), size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()) {

      for (int i = 0; i != Nblocks(); ++i)
        data_[i] = std::make_shared<Matrix>(shells_[0]->nbasis(), shells_[1]->nbasis(), true);
    }

    void compute() { compute<void*>(nullptr); }

    template<typename Value>
    void compute(const Value&) {
      static_assert(std::is_same<Value, void*>::value, "SmallInts1e::compute called illegally");
      const int s0size = shells_[0]->nbasis();
      const int s1size = shells_[1]->nbasis();
      const int a0size_inc = shells_[0]->aux_increment()->nbasis();
      const int a1size_inc = shells_[1]->aux_increment()->nbasis();
      const int a0size_dec = shells_[0]->aux_decrement() ? shells_[0]->aux_decrement()->nbasis() : 0;
      const int a1size_dec = shells_[1]->aux_decrement() ? shells_[1]->aux_decrement()->nbasis() : 0;
      const int a0 = a0size_inc + a0size_dec;
      const int a1 = a1size_inc + a1size_dec;

      constexpr int N = Batch::Nblocks();

      std::array<std::shared_ptr<Matrix>,N> unc;
      for (int n = 0; n != N; ++n)
        unc[n] = std::make_shared<Matrix>(a0, a1, true);

      {
        auto uncc = std::make_shared<Batch>(std::array<std::shared_ptr<const Shell>,2>{{shells_[0]->aux_increment(), shells_[1]->aux_increment()}}, mol_);
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(0, 0, a0size_inc, a1size_inc, uncc->data(n));
      }
      if (shells_[0]->aux_decrement() && shells_[1]->aux_decrement()) {
        auto uncc = std::make_shared<Batch>(std::array<std::shared_ptr<const Shell>,2>{{shells_[0]->aux_decrement(), shells_[1]->aux_decrement()}}, mol_);
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(a0size_inc, a1size_inc, a0size_dec, a1size_dec, uncc->data(n));
      }
      if (shells_[0]->aux_decrement()) {
        auto uncc = std::make_shared<Batch>(std::array<std::shared_ptr<const Shell>,2>{{shells_[0]->aux_decrement(), shells_[1]->aux_increment()}}, mol_);
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(a0size_inc, 0, a0size_dec, a1size_inc, uncc->data(n));
      }
      if (shells_[1]->aux_decrement()) {
        auto uncc = std::make_shared<Batch>(std::array<std::shared_ptr<const Shell>,2>{{shells_[0]->aux_increment(), shells_[1]->aux_decrement()}}, mol_);
        uncc->compute();
        for (int n = 0; n != N; ++n)
          unc[n]->copy_block(0, a1size_inc, a0size_inc, a1size_dec, uncc->data(n));
      }

      transform(unc);
    }

    std::shared_ptr<Matrix> operator[](const int i) { return data_[i]; }

    size_t size_block() const { return size_block_; }
    constexpr static int Nblocks() { return Batch::Nblocks() * 4; }

};


// in the case of finite nucleus
template<>
template<typename Value>
#ifdef LIBINT_INTERFACE
void SmallInts1e<Libint>::compute(const Value& nshells) {
#else
void SmallInts1e<ERIBatch>::compute(const Value& nshells) {
#endif
  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int a0size_inc = shells_[0]->aux_increment()->nbasis();
  const int a1size_inc = shells_[1]->aux_increment()->nbasis();
  const int a0size_dec = shells_[0]->aux_decrement() ? shells_[0]->aux_decrement()->nbasis() : 0;
  const int a1size_dec = shells_[1]->aux_decrement() ? shells_[1]->aux_decrement()->nbasis() : 0;
  const int a0 = a0size_inc + a0size_dec;
  const int a1 = a1size_inc + a1size_dec;

  auto dummy = std::make_shared<const Shell>(shells_[0]->spherical());

#ifdef LIBINT_INTERFACE
  constexpr int N = Libint::Nblocks();
#else
  constexpr int N = ERIBatch::Nblocks();
#endif
  static_assert(N == 1, "ERIBatch::Nblocks() should be 1");

  std::array<std::shared_ptr<Matrix>,N> unc;
  for (int n = 0; n != N; ++n)
    unc[n] = std::make_shared<Matrix>(a0, a1, true);

  for (auto& nshell : nshells) {
    {
#ifdef LIBINT_INTERFACE
      auto uncc = std::make_shared<Libint>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_increment(), shells_[1]->aux_increment()}});
#else
      auto uncc = std::make_shared<ERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_increment(), shells_[1]->aux_increment()}}, 2.0);
#endif
      uncc->compute();
      for (int n = 0; n != N; ++n)
        unc[n]->copy_block(0, 0, a0size_inc, a1size_inc, uncc->data(n));
    }
    if (shells_[0]->aux_decrement() && shells_[1]->aux_decrement()) {
#ifdef LIBINT_INTERFACE
      auto uncc = std::make_shared<Libint>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_decrement(), shells_[1]->aux_decrement()}});
#else
      auto uncc = std::make_shared<ERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_decrement(), shells_[1]->aux_decrement()}}, 2.0);
#endif
      uncc->compute();
      for (int n = 0; n != N; ++n)
        unc[n]->copy_block(a0size_inc, a1size_inc, a0size_dec, a1size_dec, uncc->data(n));
    }
    if (shells_[0]->aux_decrement()) {
#ifdef LIBINT_INTERFACE
      auto uncc = std::make_shared<Libint>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_decrement(), shells_[1]->aux_increment()}});
#else
      auto uncc = std::make_shared<ERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_decrement(), shells_[1]->aux_increment()}}, 2.0);
#endif
      uncc->compute();
      for (int n = 0; n != N; ++n)
        unc[n]->copy_block(a0size_inc, 0, a0size_dec, a1size_inc, uncc->data(n));
    }
    if (shells_[1]->aux_decrement()) {
#ifdef LIBINT_INTERFACE
      auto uncc = std::make_shared<Libint>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_increment(), shells_[1]->aux_decrement()}});
#else
      auto uncc = std::make_shared<ERIBatch>(std::array<std::shared_ptr<const Shell>,4>{{dummy, nshell, shells_[0]->aux_increment(), shells_[1]->aux_decrement()}}, 2.0);
#endif
      uncc->compute();
      for (int n = 0; n != N; ++n)
        unc[n]->copy_block(0, a1size_inc, a0size_inc, a1size_dec, uncc->data(n));
    }
  }

  transform(unc);
}

}

#endif
