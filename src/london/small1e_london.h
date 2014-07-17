//
// BAGEL - Parallel electron correlation program.
// Filename: small1e_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_LONDON_SMALL1E_LONDON_H
#define __SRC_LONDON_SMALL1E_LONDON_H

#include <src/math/zmatrix.h>
#include <src/integral/smallints1e_london.h>
#include <src/molecule/matrix1earray.h>

namespace bagel {

template <typename Batch>
class Small1e_London : public Matrix1eArray<4*Batch::Nblocks(), ZMatrix> {
  protected:
    void init(std::shared_ptr<const Molecule> mol) override {
      std::list<std::shared_ptr<const Shell>> shells;
      for (auto& i : mol->atoms())
        shells.insert(shells.end(), i->shells().begin(), i->shells().end());

      // TODO thread, parallel
      int o0 = 0;
      for (auto& a0 : shells) {
        int o1 = 0;
        for (auto& a1 : shells) {
          std::array<std::shared_ptr<const Shell>,2> input = {{a1, a0}};
          computebatch(input, o0, o1, mol);
          o1 += a1->nbasis();
        }
        o0 += a0->nbasis();
      }
    }


  public:
    Small1e_London(const std::shared_ptr<const Molecule> mol) : Matrix1eArray<4*Batch::Nblocks(), ZMatrix>(mol) {
      init(mol);
    }

    void computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, std::shared_ptr<const Molecule> mol) {
      // input = [b1, b0]
      assert(input.size() == 2);
      const int dimb1 = input[0]->nbasis();
      const int dimb0 = input[1]->nbasis();
      SmallInts1e_London<Batch> batch(input, mol);
      batch.compute();

      for (int i = 0; i != this->Nblocks(); ++i)
        this->matrices_[i]->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch[i]);
    }

    void print(const std::string name = "") const override {
      Matrix1eArray<4*Batch::Nblocks(), ZMatrix>::print(name.empty() ? "Small1e" : name);
    }
};

template<> void Small1e_London<ComplexERIBatch>::computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, std::shared_ptr<const Molecule>);

}

#endif
