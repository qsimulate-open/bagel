//
// BAGEL - Parallel electron correlation program.
// Filename: df_london.h
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

#ifndef __SRC_DF_DF_LONDON_H
#define __SRC_DF_DF_LONDON_H

#include <src/molecule/atom.h>
#include <src/parallel/resources.h>
#include <src/math/zmatrix.h>
#include <src/df/paralleldf_london.h>
#include <src/df/dfblock_london.h>
#include <src/df/dfinttask.h>
#include <src/df/dfinttask_old.h>
#include <src/integral/comprys/complexeribatch.h>
#include <src/util/timer.h>
#include <src/util/taskqueue.h>

namespace bagel {

class DFHalfDist_London;
class DFFullDist_London;


class DFDist_London : public ParallelDF_London {
  friend class DFIntTask_OLD<DFDist_London>;
  protected:
    std::pair<const double*, std::shared_ptr<RysIntegral<double, Int_t::Standard>>> compute_batch(std::array<std::shared_ptr<const Shell>,4>& input);

    std::shared_ptr<const StaticDist> make_table(const size_t nmax);
    // compute 2-index integrals ERI
    void compute_2index(const std::vector<std::shared_ptr<const Shell>>&, const double thresh, const bool compute_inv);

    std::tuple<int, std::vector<std::shared_ptr<const Shell>>> get_ashell(const std::vector<std::shared_ptr<const Shell>>& all);

  public:
    DFDist_London(const int nbas, const int naux, const std::shared_ptr<DFBlock_London> block = nullptr, std::shared_ptr<const ParallelDF_London> df = nullptr, std::shared_ptr<Matrix> data2 = nullptr)
      : ParallelDF_London(naux, nbas, nbas, df, data2) {
      if (block)
        block_.push_back(block);
    }

    DFDist_London(const std::shared_ptr<const ParallelDF_London> df) : ParallelDF_London(df->naux(), df->nindex1(), df->nindex2(), df) { }

    bool has_2index() const { return data2_.get() != nullptr; }
    size_t nbasis0() const { return nindex2_; }
    size_t nbasis1() const { return nindex1_; }
    size_t naux() const { return naux_; }

    void add_direct_product(std::shared_ptr<const ZMatrix> a, std::shared_ptr<const ZMatrix> b, const double fac)
       { add_direct_product(std::vector<std::shared_ptr<const ZMatrix>>{a}, std::vector<std::shared_ptr<const ZMatrix>>{b}, fac); }
    void add_direct_product(std::vector<std::shared_ptr<const ZMatrix>> a, std::vector<std::shared_ptr<const ZMatrix>> b, const double fac);

    // compute half transforms; c is dimensioned by nbasis_;
    std::shared_ptr<DFHalfDist_London> compute_half_transform(const ZMatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<DFHalfDist_London> compute_half_transform(std::shared_ptr<T> c) const { return compute_half_transform(*c); }

    // compute half transform using the third index. You get DFHalfDist with gamma/i/s (i.e., index are reordered)
    std::shared_ptr<DFHalfDist_London> compute_half_transform_swap(const ZMatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<DFHalfDist_London> compute_half_transform_swap(std::shared_ptr<T> c) const { return compute_half_transform_swap(*c); }

    std::shared_ptr<DFDist_London> copy() const;
    std::shared_ptr<DFDist_London> clone() const;

    // split up smalleri integrals into 6 dfdist objects
    std::vector<std::shared_ptr<const DFDist_London>> split_blocks() const {
      std::vector<std::shared_ptr<const DFDist_London>> out;
      assert(nindex1_ == nindex2_);
      for (auto& i : block_)
        out.push_back(std::make_shared<const DFDist_London>(nindex1_, naux_, i, df_, data2_));
      return out;
    }
};


template<class TBatch>
class DFDist_London_ints : public DFDist_London {
  protected:
    void compute_3index(const std::vector<std::shared_ptr<const Shell>>& ashell,
                        const std::vector<std::shared_ptr<const Shell>>& b1shell,
                        const std::vector<std::shared_ptr<const Shell>>& b2shell,
                        const size_t asize, const size_t b1size, const size_t b2size,
                        const size_t astart, const double thresh, const bool compute_inv) {
      Timer time;

      // making a task list
      TaskQueue<DFIntTask<TBatch, TBatch::Nblocks(), DFBlock_London, std::complex<double>>> tasks(b1shell.size()*b2shell.size()*ashell.size());

      auto i3 = std::make_shared<const Shell>(ashell.front()->spherical());

      // due to performance issue, we need to reshape it to array
      std::array<std::shared_ptr<DFBlock_London>,TBatch::Nblocks()> blk;
      for (int i = 0; i != TBatch::Nblocks(); ++i) blk[i] = block_[i];

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
    DFDist_London_ints(const int nbas, const int naux, const std::vector<std::shared_ptr<const Atom>>& atoms, const std::vector<std::shared_ptr<const Atom>>& aux_atoms,
                const double thr, const bool inverse, const double dum, const bool average = false) : DFDist_London(nbas, naux) {

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
      for (int i = 0; i != TBatch::Nblocks(); ++i)
        block_.push_back(std::make_shared<DFBlock_London>(adist_shell, adist_averaged, asize, b1size, b2size, astart, 0, 0));

      // 3-index integrals
      compute_3index(myashell, b1shell, b2shell, asize, b1size, b2size, astart, thr, inverse);
      // 2-index integrals
      compute_2index(ashell, thr, inverse);
      // 3-index integrals, post process
      if (average)
        average_3index();
    }

};


class DFHalfDist_London : public ParallelDF_London {
  protected:

  public:
    DFHalfDist_London(const std::shared_ptr<const ParallelDF_London> df, const int nocc) : ParallelDF_London(df->naux(), nocc, df->nindex2(), df) { }

    size_t nocc() const { return nindex1_; }
    size_t nbasis() const { return nindex2_; }

    std::shared_ptr<DFFullDist_London> compute_second_transform(const ZMatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<DFFullDist_London> compute_second_transform(std::shared_ptr<T> c) const { return compute_second_transform(*c); }

    std::shared_ptr<DFDist_London> back_transform(const ZMatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<DFDist_London> back_transform(std::shared_ptr<T> c) const { return back_transform(*c); }

    std::shared_ptr<DFHalfDist_London> copy() const;
    std::shared_ptr<DFHalfDist_London> clone() const;

    void rotate_occ(const std::shared_ptr<const ZMatrix> d);
    std::shared_ptr<DFHalfDist_London> apply_density(const std::shared_ptr<const ZMatrix> d) const;

    std::shared_ptr<ZMatrix> compute_Kop_1occ(const std::shared_ptr<const ZMatrix> den, const double a) const;

    std::shared_ptr<DFHalfDist_London> apply_J() const { return apply_J(df_->data2()); }
    std::shared_ptr<DFHalfDist_London> apply_JJ() const { return apply_J(std::make_shared<ZMatrix>(*df_->data2()**df_->data2())); }
    std::shared_ptr<DFHalfDist_London> apply_J(const std::shared_ptr<const DFDist_London> d) const { return apply_J(d->data2()); }
    std::shared_ptr<DFHalfDist_London> apply_JJ(const std::shared_ptr<const DFDist_London> d) const { return apply_J(std::make_shared<ZMatrix>(*d->data2()**d->data2())); }
    std::shared_ptr<DFHalfDist_London> apply_J(const std::shared_ptr<const ZMatrix> o) const;

};


class DFFullDist_London : public ParallelDF_London {
  protected:
    std::shared_ptr<DFFullDist_London> apply_J(const std::shared_ptr<const ZMatrix> o) const;

  public:
    DFFullDist_London(const std::shared_ptr<const ParallelDF_London> df, const int nocc1, const int nocc2) : ParallelDF_London(df->naux(), nocc1, nocc2, df) { }

    int nocc1() const { return nindex1_; }
    int nocc2() const { return nindex2_; }

    std::shared_ptr<DFFullDist_London> copy() const;
    std::shared_ptr<DFFullDist_London> clone() const;

    std::shared_ptr<DFHalfDist_London> back_transform(const ZMatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<DFHalfDist_London> back_transform(std::shared_ptr<T> c) const { return back_transform(*c); }

    void rotate_occ1(const std::shared_ptr<const ZMatrix> d);

/*
    // 2RDM contractions
    // special function for RHF
    std::shared_ptr<DFFullDist_London> apply_closed_2RDM(const double scale_ex = 1.0) const;
    // special function for UHF
    std::shared_ptr<DFFullDist_London> apply_uhf_2RDM(const double* rdma, const double* rdmb) const;
    // general case with closed orbitals
    std::shared_ptr<DFFullDist_London> apply_2rdm(const double* rdm, const double* rdm1, const int nclosed, const int nact) const;
    // general case without closed orbitals
    std::shared_ptr<DFFullDist_London> apply_2rdm(const double* rdm) const;
*/

    std::shared_ptr<ZMatrix> form_4index_1fixed(const std::shared_ptr<const DFFullDist_London> o, const double a, const size_t n) const;

    // utility functions
    std::shared_ptr<ZMatrix> form_aux_2index_apply_J(const std::shared_ptr<const DFFullDist_London> o, const double a) const;
    void add_product(const std::shared_ptr<const DFFullDist_London>, const std::shared_ptr<const ZMatrix>, const int jdim, const size_t offset, const double fac = 1.0);

    std::shared_ptr<DFFullDist_London> apply_J() const { return apply_J(df_->data2()); }
    std::shared_ptr<DFFullDist_London> apply_JJ() const { return apply_J(std::make_shared<ZMatrix>(*df_->data2()**df_->data2())); }
    std::shared_ptr<DFFullDist_London> apply_J(const std::shared_ptr<const ParallelDF_London> d) const { return apply_J(d->data2()); }
    std::shared_ptr<DFFullDist_London> apply_JJ(const std::shared_ptr<const ParallelDF_London> d) const { return apply_J(std::make_shared<ZMatrix>(*d->data2()**d->data2())); }

    std::shared_ptr<DFFullDist_London> swap() const;
};

}


#endif
