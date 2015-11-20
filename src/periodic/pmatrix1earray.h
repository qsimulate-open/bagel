//
// BAGEL - Parallel electron correlation program.
// Filename: pmatrix1earray.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_PERIODIC_PMATRIX1EARRAY_H
#define __SRC_PERIODIC_PMATRIX1EARRAY_H

#include <src/util/constants.h>
#include <src/periodic/pmatrix1e.h>


namespace bagel {

template <int N>
class PMatrix1eArrayTask;

// N component periodic 1e integrals
template <int N>
class PMatrix1eArray {
  friend class PMatrix1eArrayTask<N>;
  protected:
    std::array<std::shared_ptr<PData>, N> pdata_blocks_;

    virtual void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int, std::shared_ptr<const Lattice>, const int) = 0;
    virtual void init(std::shared_ptr<const Lattice>);

  public:
    PMatrix1eArray() { }
    PMatrix1eArray(const std::shared_ptr<const Lattice>);
    virtual ~PMatrix1eArray() { }

    std::shared_ptr<PData>& data(const int i) { return pdata_blocks_[i]; }
    std::shared_ptr<const PData> data(const int i) const { return pdata_blocks_[i]; }

    PData& operator[](const int i) { return *pdata_blocks_[i]; }
    const PData& operator[](const int i) const { return *pdata_blocks_[i]; }

    constexpr static int Nblocks() { return N; }

    void fill_upper_conjg() { for (int i = 0 ; i < N; ++i) pdata_blocks_[i]->fill_upper_conjg(); }
};



template <int N>
PMatrix1eArray<N>::PMatrix1eArray(const std::shared_ptr<const Lattice> lattice) {
  static_assert(N > 0, "PMatrix1eArray should be constructed with N > 0");
  for(int i = 0; i < N; ++i) {
    pdata_blocks_[i] = std::make_shared<PData>(lattice->primitive_cell()->nbasis(), lattice->num_lattice_vectors());
  }
}


template <int N>
void PMatrix1eArray<N>::init(std::shared_ptr<const Lattice> lattice) {

  std::shared_ptr<const Geometry> cell0 = lattice->primitive_cell();
  const size_t nshell = accumulate(cell0->atoms().begin(), cell0->atoms().end(), 0, [](int r, std::shared_ptr<const Atom> p) { return r+p->nshell(); });
  TaskQueue<PMatrix1eArrayTask<N>> task(nshell * nshell * lattice->extent());

  int g = 0;
  int u = 0;
  /* loop over all cells in direct space */
  for (auto& disp : lattice->lattice_vectors()) {
    auto cell = std::make_shared<const Geometry>(*(lattice->primitive_cell()), disp);

    size_t oa0 = 0;
    for (auto a0 = cell0->atoms().begin(); a0 != cell0->atoms().end(); ++a0) { /* cell 0 */
      size_t oa1 = 0;
      for (auto a1 = cell->atoms().begin(); a1 != cell->atoms().end(); ++a1) { /* cell g */

        size_t ob0 = oa0;
        for (auto& b0 : (*a0)->shells()) {
          size_t ob1 = oa1;
          for (auto& b1 : (*a1)->shells()) {
            if (u++ % mpi__->size() == mpi__->rank())
              task.emplace_back(std::array<std::shared_ptr<const Shell>,2>{{b1, b0}}, ob0, ob1, lattice, this, g);
            ob1 += b1->nbasis();
          }
          ob0 += b0->nbasis();
        }
        oa1 += (*a1)->nbasis();
      }
      oa0 += (*a0)->nbasis();
    }
    ++g;
  }

  task.compute();
  for (auto& block : pdata_blocks_) block->allreduce();
}


template <int N>
class PMatrix1eArrayTask {
  protected:
    PMatrix1eArray<N>* parent_;
    size_t ob0, ob1;
    std::array<std::shared_ptr<const Shell>, 2> basis_;
    std::shared_ptr<const Lattice> lattice_;
    int block_;
  public:
    PMatrix1eArrayTask<N>(std::array<std::shared_ptr<const Shell>,2> sh, size_t b0, size_t b1, std::shared_ptr<const Lattice> l,
                          PMatrix1eArray<N>* p, int blk)
      : parent_(p), ob0(b0), ob1(b1), basis_(sh), lattice_(l), block_(blk) { }
    void compute() const { parent_->computebatch(basis_, ob0, ob1, lattice_, block_); }
};

}

#endif
