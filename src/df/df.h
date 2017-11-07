//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: df.h
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

#ifndef __SRC_DF_DF_H
#define __SRC_DF_DF_H

#include <src/df/paralleldf.h>
#include <src/molecule/atom.h>

namespace bagel {

class DFHalfDist;
class DFFullDist;

class DFDist : public ParallelDF {
  friend class DFIntTask_OLD<DFDist>;
  friend class PDFIntTask_2index;
  protected:
    std::pair<const double*, std::shared_ptr<RysInt>> compute_batch(std::array<std::shared_ptr<const Shell>,4>& input);

    std::shared_ptr<const StaticDist> make_table(const size_t nmax);
    // compute 2-index integrals ERI
    void compute_2index(const std::vector<std::shared_ptr<const Shell>>&, const double thresh, const bool compute_inv);

    std::tuple<int, std::vector<std::shared_ptr<const Shell>>> get_ashell(const std::vector<std::shared_ptr<const Shell>>& all);

  public:
    DFDist(const int nbas, const int naux, const std::shared_ptr<DFBlock> block = nullptr, std::shared_ptr<const ParallelDF> df = nullptr, std::shared_ptr<Matrix> data2 = nullptr,
           const bool serial = false) : ParallelDF(naux, nbas, nbas, df, data2, serial) {
      if (block)
        block_.push_back(block);
    }

    DFDist(const std::shared_ptr<const ParallelDF> df) : ParallelDF(df->naux(), df->nindex1(), df->nindex2(), df) { }

    bool has_2index() const { return data2_.get() != nullptr; }
    size_t nbasis0() const { return nindex2_; }
    size_t nbasis1() const { return nindex1_; }
    size_t naux() const { return naux_; }

    // compute half transforms; c is dimensioned by nbasis_;
    std::shared_ptr<DFHalfDist> compute_half_transform(const MatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<DFHalfDist> compute_half_transform(std::shared_ptr<T> c) const { return compute_half_transform(*c); }

    // compute half transform using the third index. You get DFHalfDist with gamma/i/s (i.e., index are reordered)
    std::shared_ptr<DFHalfDist> compute_half_transform_swap(const MatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<DFHalfDist> compute_half_transform_swap(std::shared_ptr<T> c) const { return compute_half_transform_swap(*c); }

    std::shared_ptr<DFDist> copy() const;
    std::shared_ptr<DFDist> clone() const;

    // split up smalleri integrals into 6 dfdist objects
    virtual std::vector<std::shared_ptr<const DFDist>> split_blocks() const {
      std::vector<std::shared_ptr<const DFDist>> out;
      assert(nindex1_ == nindex2_);
      for (auto& i : block_)
        out.push_back(std::make_shared<const DFDist>(nindex1_, naux_, i, df_, data2_));
      return out;
    }
};


template<class TBatch>
class DFDist_ints : public DFDist {
  protected:
    void compute_3index(const std::vector<std::shared_ptr<const Shell>>& ashell,
                        const std::vector<std::shared_ptr<const Shell>>& b1shell,
                        const std::vector<std::shared_ptr<const Shell>>& b2shell,
                        const size_t asize, const size_t b1size, const size_t b2size,
                        const size_t astart, const double thresh, const bool compute_inv) {
      Timer time;

      // making a task list
      TaskQueue<DFIntTask<TBatch,TBatch::Nblocks()>> tasks(b1shell.size()*b2shell.size()*ashell.size());

      auto i3 = std::make_shared<const Shell>(ashell.front()->spherical());

      // due to performance issue, we need to reshape it to array
      std::array<std::shared_ptr<DFBlock>,TBatch::Nblocks()> blk;
      for (int i = 0; i != TBatch::Nblocks(); ++i) blk[i] = block_[i];

      int j2 = 0;
      for (auto& i2 : b2shell) {
        int j1 = 0;
        for (auto& i1 : b1shell) {
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
    DFDist_ints(const int nbas, const int naux, const std::vector<std::shared_ptr<const Atom>>& atoms, const std::vector<std::shared_ptr<const Atom>>& aux_atoms,
                const double thr, const bool inverse, const double dum, const bool average = false, const std::shared_ptr<Matrix> data2 = nullptr, const bool serial = false)
      : DFDist(nbas, naux, nullptr, nullptr, nullptr, serial) {

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
        block_.push_back(std::make_shared<DFBlock>(adist_shell, adist_averaged, asize, b1size, b2size, astart, 0, 0));

      // 3-index integrals
      compute_3index(myashell, b1shell, b2shell, asize, b1size, b2size, astart, thr, inverse);

      // 2-index integrals
      if (data2)
        data2_ = data2;
      else
        compute_2index(ashell, thr, inverse);

      // 3-index integrals, post process
      if (average)
        average_3index();
    }

};


class DFHalfDist : public ParallelDF {
  protected:

  public:
    DFHalfDist(const std::shared_ptr<const ParallelDF> df, const int nocc) : ParallelDF(df->naux(), nocc, df->nindex2(), df) { }

    size_t nocc() const { return nindex1_; }
    size_t nbasis() const { return nindex2_; }

    std::shared_ptr<DFFullDist> compute_second_transform(const MatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<DFFullDist> compute_second_transform(std::shared_ptr<T> c) const { return compute_second_transform(*c); }

    std::shared_ptr<DFDist> back_transform(const MatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<DFDist> back_transform(std::shared_ptr<T> c) const { return back_transform(*c); }

    std::shared_ptr<DFHalfDist> copy() const;
    std::shared_ptr<DFHalfDist> clone() const;

    std::shared_ptr<DFHalfDist> merge_b1(std::shared_ptr<DFHalfDist> o) const;
    std::shared_ptr<DFHalfDist> slice_b1(const int slice_start, const int slice_size) const;

    std::shared_ptr<DFHalfDist> transform_occ(const std::shared_ptr<const Matrix> d) const;
    std::shared_ptr<DFHalfDist> apply_density(const std::shared_ptr<const Matrix> d) const;

    std::shared_ptr<Matrix> compute_Kop_1occ(const std::shared_ptr<const Matrix> den, const double a) const;

    std::shared_ptr<DFHalfDist> apply_J() const { return apply_J(df_->data2()); }
    std::shared_ptr<DFHalfDist> apply_JJ() const { return apply_J(std::make_shared<Matrix>(*df_->data2()**df_->data2())); }
    std::shared_ptr<DFHalfDist> apply_J(const std::shared_ptr<const DFDist> d) const { return apply_J(d->data2()); }
    std::shared_ptr<DFHalfDist> apply_JJ(const std::shared_ptr<const DFDist> d) const { return apply_J(std::make_shared<Matrix>(*d->data2()**d->data2())); }
    std::shared_ptr<DFHalfDist> apply_J(const std::shared_ptr<const Matrix> o) const;

};


class DFFullDist : public ParallelDF {
  protected:
    std::shared_ptr<DFFullDist> apply_J(const std::shared_ptr<const Matrix> o) const;

  public:
    DFFullDist(const std::shared_ptr<const ParallelDF> df, const int nocc1, const int nocc2) : ParallelDF(df->naux(), nocc1, nocc2, df) { }

    int nocc1() const { return nindex1_; }
    int nocc2() const { return nindex2_; }

    std::shared_ptr<DFFullDist> copy() const;
    std::shared_ptr<DFFullDist> clone() const;

    std::shared_ptr<DFHalfDist> back_transform(const MatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<DFHalfDist> back_transform(std::shared_ptr<T> c) const { return back_transform(*c); }

    std::shared_ptr<DFFullDist> transform_occ1(const std::shared_ptr<const Matrix> d) const;

    // 2RDM contractions
    // special function for RHF
    std::shared_ptr<DFFullDist> apply_closed_2RDM(const double scale_ex = 1.0) const;
    // special function for UHF
    std::shared_ptr<DFFullDist> apply_uhf_2RDM(const btas::Tensor2<double>& rdma, const btas::Tensor2<double>& rdmb) const;
    // general case with closed orbitals
    std::shared_ptr<DFFullDist> apply_2rdm(const btas::Tensor4<double>& rdm, const btas::Tensor2<double>& rdm1, const int nclosed, const int nact) const;
    // general case without closed orbitals
    std::shared_ptr<DFFullDist> apply_2rdm(const btas::Tensor4<double>& rdm) const;
    // general case with closed orbitals, transition case
    std::shared_ptr<DFFullDist> apply_2rdm_tran(const btas::Tensor4<double>& rdm, const btas::Tensor2<double>& rdm1, const int nclosed, const int nact) const;

    // returns the 4-index integrals with fixed index n
    std::shared_ptr<Matrix> form_4index_1fixed(const std::shared_ptr<const DFFullDist> o, const double a, const size_t n) const;
    // returns the diagonal part of the 4-index integrals [o1v0] = (o1v0|g)(g|o1v0)
    std::shared_ptr<Matrix> form_4index_diagonal() const;
    // returns the diagonal part of the 4-index integrals [o1o2v0] = (o1v0|g)(g|o2v0)
    std::shared_ptr<Matrix> form_4index_diagonal_part() const;

    // utility functions
    std::shared_ptr<Matrix> form_aux_2index_apply_J(const std::shared_ptr<const DFFullDist> o, const double a) const;

    std::shared_ptr<DFFullDist> apply_J() const { return apply_J(df_->data2()); }
    std::shared_ptr<DFFullDist> apply_JJ() const { return apply_J(std::make_shared<Matrix>(*df_->data2()**df_->data2())); }
    std::shared_ptr<DFFullDist> apply_J(const std::shared_ptr<const ParallelDF> d) const { return apply_J(d->data2()); }
    std::shared_ptr<DFFullDist> apply_JJ(const std::shared_ptr<const ParallelDF> d) const { return apply_J(std::make_shared<Matrix>(*d->data2()**d->data2())); }

    std::shared_ptr<DFFullDist> swap() const;
};

}


#endif
