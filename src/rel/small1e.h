//
// BAGEL - Parallel electron correlation program.
// Filename: small1e.h
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


#ifndef __SRC_REL_SMALL1E_H
#define __SRC_REL_SMALL1E_H

#include <src/math/zmatrix.h>
#include <src/integral/smallints1e.h>
#include <src/molecule/matrix1earray.h>

namespace bagel {

template <typename Batch>
class Small1e : public Matrix1eArray<4*Batch::Nblocks()> {
  protected:
    void init(std::shared_ptr<const Molecule> mol) override {
      std::list<std::shared_ptr<const Shell>> shells;
      for (auto& i : mol->atoms())
        shells.insert(shells.end(), i->shells().begin(), i->shells().end());

      const size_t nshell = accumulate(mol->atoms().begin(), mol->atoms().end(), 0, [](int r, std::shared_ptr<const Atom> p) { return r+p->nshell(); });
      TaskQueue<Matrix1eArrayTask<4*Batch::Nblocks()>> task(nshell*(nshell+1)/2);

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

    void computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, std::shared_ptr<const Molecule> mol) override {
      // input = [b1, b0]
      assert(input.size() == 2);
      const int dimb1 = input[0]->nbasis();
      const int dimb0 = input[1]->nbasis();
      SmallInts1e<Batch> batch(input, mol);
      batch.compute();

      for (int i = 0; i != this->Nblocks(); ++i)
        this->matrices_[i]->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch[i]);
    }

  public:
    Small1e(const std::shared_ptr<const Molecule> mol) : Matrix1eArray<4*Batch::Nblocks()>(mol) {
      init(mol);
    }

    void print(const std::string name = "") const override {
      Matrix1eArray<4*Batch::Nblocks()>::print(name.empty() ? "Small1e" : name);
    }
};

template<> void Small1e<ERIBatch>::computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, std::shared_ptr<const Molecule>);

}

#endif
