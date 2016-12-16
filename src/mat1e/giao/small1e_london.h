//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: small1e_london.h
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


#ifndef __SRC_LONDON_SMALL1E_LONDON_H
#define __SRC_LONDON_SMALL1E_LONDON_H

#include <src/util/math/zmatrix.h>
#include <src/integral/smallints1e_london.h>
#include <src/integral/comprys/complexnaibatch.h>
#include <src/mat1e/matrix1earray.h>
#include <src/util/unpack.h>

namespace bagel {

// "Args" will be saved and passed to each batch of integrals
// If the shared_ptr<const Molecule> is needed, it can be redundantly supplied in Args or just passed directly in a specialized computebatch()
template<typename Batch, typename ...Args>
class Small1e_London : public Matrix1eArray<4*Batch::Nblocks(), ZMatrix> {
  protected:
    const std::tuple<Args...> args_;

    void init(std::shared_ptr<const Molecule> mol) override {
      std::list<std::shared_ptr<const Shell>> shells;
      for (auto& i : mol->atoms())
        shells.insert(shells.end(), i->shells().begin(), i->shells().end());

      const size_t nshell = accumulate(mol->atoms().begin(), mol->atoms().end(), 0, [](int r, std::shared_ptr<const Atom> p) { return r+p->nshell(); });
      TaskQueue<Matrix1eArrayTask<4*Batch::Nblocks(), ZMatrix>> task(nshell*(nshell+1)/2);

      int o0 = 0;
      int u = 0;
      for (auto& a0 : shells) {
        int o1 = 0;
        for (auto& a1 : shells) {
          if (u++ % mpi__->size() == mpi__->rank()) {
            std::array<std::shared_ptr<const Shell>,2> input = {{a1, a0}};
            task.emplace_back(input, o0, o1, mol, this);
          }
          o1 += a1->nbasis();
        }
        o0 += a0->nbasis();
      }
      task.compute();
      for (auto& i : this->matrices_) i->allreduce();
    }

    // Unpack variadic template arguments here
    template <int ...S>
    SmallInts1e_London<Batch, Args...> get_batch(std::array<std::shared_ptr<const Shell>,2> input, seq<S...>) {
      SmallInts1e_London<Batch, Args...> out(input, std::get<S>(args_) ...);
      return out;
    }

    void computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, std::shared_ptr<const Molecule> mol) override {
      // input = [b1, b0]
      assert(input.size() == 2);
      const int dimb1 = input[0]->nbasis();
      const int dimb0 = input[1]->nbasis();
      SmallInts1e_London<Batch, Args...> batch = get_batch(input, typename gens<sizeof...(Args)>::type());
      batch.compute();

      for (int i = 0; i != this->Nblocks(); ++i)
        this->matrices_[i]->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch[i]);
    }

  public:
    Small1e_London(const std::shared_ptr<const Molecule> mol, Args... args) : Matrix1eArray<4*Batch::Nblocks(), ZMatrix>(mol), args_(args...) {
      init(mol);
    }

    void print(const std::string name = "", const int len = 10) const override {
      Matrix1eArray<4*Batch::Nblocks(), ZMatrix>::print(name.empty() ? "Small1e" : name, len);
    }
};

template<> void Small1e_London<ComplexNAIBatch>::computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, std::shared_ptr<const Molecule>);
template<> void Small1e_London<ComplexERIBatch>::computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, std::shared_ptr<const Molecule>);

}

#endif
