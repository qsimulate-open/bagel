//
// BAGEL - Parallel electron correlation program.
// Filename: complexdf.h
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


#ifndef __SRC_DF_COMPLEXDF_H
#define __SRC_DF_COMPLEXDF_H

#include <src/df/df.h>
#include <src/df/complexdfinttask.h>
#include <src/df/complexparalleldf.h>
#include <src/math/zmatrix.h>

namespace bagel {

class ComplexDFHalfDist;

class ComplexDFDist : public ComplexParallelDF {
  protected:
    std::shared_ptr<const StaticDist> make_table(const size_t nmax);
    std::tuple<int, std::vector<std::shared_ptr<const Shell>>> get_ashell(const std::vector<std::shared_ptr<const Shell>>& all);
    void compute_2index(const std::vector<std::shared_ptr<const Shell>>&, const double thresh, const bool compute_inv);

  public:
    ComplexDFDist(const int nbas, const int naux, const int nblock, const std::shared_ptr<DFBlock> block = nullptr,
                  std::shared_ptr<const ComplexParallelDF> df = nullptr, std::shared_ptr<Matrix> data2 = nullptr);

    ComplexDFDist(const std::shared_ptr<const ComplexParallelDF> df);
    ComplexDFDist(const std::array<std::shared_ptr<DFBlock>,2> block, const std::shared_ptr<const ComplexParallelDF> df, const std::shared_ptr<Matrix> d2);

    bool has_2index() const { return data2_.get() != nullptr; }
    size_t nbasis0() const { return nindex2_; }
    size_t nbasis1() const { return nindex1_; }
    size_t naux() const { return naux_; }

    // compute half transforms; c is dimensioned by nbasis_;
    std::shared_ptr<ComplexDFHalfDist> compute_half_transform(const ZMatView c) const;

    // compute half transform using the third index. You get DFHalfDist with gamma/i/s (i.e., index are reordered)
    std::shared_ptr<ComplexDFHalfDist> compute_half_transform_swap(const ZMatView c) const;

    // split up smalleri integrals into 6 dfdist objects
    std::vector<std::shared_ptr<const ComplexDFDist>> split_blocks() const {
      assert(nindex1_ == nindex2_);
      std::vector<std::shared_ptr<DFBlock>> outr;
      std::vector<std::shared_ptr<DFBlock>> outi;
      std::vector<std::shared_ptr<const ComplexDFDist>> out;
      for (int i=0; i!=nblock_; i++) {
        outr.push_back(dfdata_[0]->block_[i]);
        outi.push_back(dfdata_[1]->block_[i]);
        std::array<std::shared_ptr<DFBlock>,2> in = {{ outr[i], outi[i] }};
        out.push_back(std::make_shared<const ComplexDFDist>(in, df_ ? df_ : shared_from_this() , data2_));
      }
      return out;
    }
};

template<class TBatch>
class ComplexDFDist_ints : public ComplexDFDist {
  protected:
    void compute_3index(const std::vector<std::shared_ptr<const Shell>>& ashell,
                        const std::vector<std::shared_ptr<const Shell>>& b1shell,
                        const std::vector<std::shared_ptr<const Shell>>& b2shell,
                        const size_t asize, const size_t b1size, const size_t b2size,
                        const size_t astart, const double thresh, const bool compute_inv) {
      Timer time;

      // making a task list
      TaskQueue<ComplexDFIntTask<TBatch,2*TBatch::Nblocks()>> tasks(b1shell.size()*b2shell.size()*ashell.size());

      auto i3 = std::make_shared<const Shell>(ashell.front()->spherical());

      // due to performance issue, we need to reshape it to array
      std::array<std::shared_ptr<DFBlock>,2*TBatch::Nblocks()> blk;
      for (int i = 0; i != TBatch::Nblocks(); ++i) {
        const int i2 = 2*i;
        blk[i2] =   get_real()->block_[i];
        blk[i2+1] = get_imag()->block_[i];
      }

      int j2 = 0;
      for (auto& i2 : b2shell) {
        int j1 = 0;
        for (auto& i1 : b1shell) {
          // TODO careful
          // Since we use a real auxiliary basis set, half the 3-index integrals will be complex conjugates of the other half
          if (TBatch::Nblocks() > 1 || j1 <= j2) {
            int j0 = 0;
            for (auto& i0 : ashell) {
              tasks.emplace_back((std::array<std::shared_ptr<const Shell>,4>{{i3, i0, i1, i2}}), (std::array<int,3>{{j2, j1, j0}}), blk);
              j0 += i0->nbasis();
            }
          }
          j1 += i1->nbasis();
        }
        j2 += i2->nbasis();
      }
      time.tick_print("3-index ints prep");
      tasks.compute();
      time.tick_print("3-index ints");

    }

  public:
    ComplexDFDist_ints(const int nbas, const int naux, const std::vector<std::shared_ptr<const Atom>>& atoms, const std::vector<std::shared_ptr<const Atom>>& aux_atoms,
                const double thr, const bool inverse, const double dum, const bool average = false) : ComplexDFDist(nbas, naux, TBatch::Nblocks()) {

      // 3index Integral is now made in DFBlock.
      std::vector<std::shared_ptr<const Shell>> ashell, b1shell, b2shell;
      for (auto& i : aux_atoms) ashell.insert(ashell.end(), i->shells().begin(), i->shells().end());
      for (auto& i : atoms)     b1shell.insert(b1shell.end(), i->shells().begin(), i->shells().end());
      for (auto& i : atoms)     b2shell.insert(b2shell.end(), i->shells().begin(), i->shells().end());

      // distribute auxiliary shells to each nodes
      int astart;
      std::vector<std::shared_ptr<const Shell>> myashell;
      std::tie(astart, myashell) = get_ashell(ashell);

      std::shared_ptr<const StaticDist> adist_shell = make_table(astart);
      std::shared_ptr<const StaticDist> adist_averaged = std::make_shared<const StaticDist>(naux_, mpi__->size());

      // make empty dfblocks
      const size_t asize  = std::accumulate(myashell.begin(),myashell.end(),0, [](const int& i, const std::shared_ptr<const Shell>& o) { return i+o->nbasis(); });
      const size_t b1size = std::accumulate(b1shell.begin(), b1shell.end(), 0, [](const int& i, const std::shared_ptr<const Shell>& o) { return i+o->nbasis(); });
      const size_t b2size = std::accumulate(b2shell.begin(), b2shell.end(), 0, [](const int& i, const std::shared_ptr<const Shell>& o) { return i+o->nbasis(); });

      for (int i = 0; i != TBatch::Nblocks(); i++) {
        dfdata_[0]->block_.push_back(std::make_shared<DFBlock>(adist_shell, adist_averaged, asize, b1size, b2size, astart, 0, 0));
        dfdata_[1]->block_.push_back(std::make_shared<DFBlock>(adist_shell, adist_averaged, asize, b1size, b2size, astart, 0, 0));
      }

      // 3-index integrals
      compute_3index(myashell, b1shell, b2shell, asize, b1size, b2size, astart, thr, inverse);
      // 2-index integrals
      compute_2index(ashell, thr, inverse);
      // 3-index integrals, post process
      if (average)
        average_3index();
    }

};


// TODO Functions not needed for RHF have not been implemented (transform_second, back_transform, etc.)
class ComplexDFHalfDist : public ComplexParallelDF {
  protected:

  public:
    //ComplexDFHalfDist(const std::shared_ptr<const ComplexParallelDF> df, const int nocc) : ComplexParallelDF(df->naux(), nocc, df->nindex2(), df->nblock(), df) { }
    ComplexDFHalfDist(std::shared_ptr<DFHalfDist> rdf, std::shared_ptr<DFHalfDist> idf, const std::shared_ptr<const ComplexParallelDF>);

    size_t nocc() const { return nindex1_; }
    size_t nbasis() const { return nindex2_; }

    std::shared_ptr<ComplexDFHalfDist> apply_J() const { return apply_J(df_->data2_real()); }
    //std::shared_ptr<ComplexDFHalfDist> apply_J(const std::shared_ptr<const ComplexDFDist> d) const { return apply_J(d->data2_real()); }
    std::shared_ptr<ComplexDFHalfDist> apply_J(const std::shared_ptr<const Matrix> o) const;

};


}

#endif
