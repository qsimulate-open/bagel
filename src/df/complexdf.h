//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexdf.h
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

// Specialized class of DF integral tensors for complex basis functions, currently only used for GIAO
// The design of this class is not entirely satisfactory, see issue #50 on github.com/nubakery/BAGEL


#ifndef __SRC_DF_COMPLEXDF_H
#define __SRC_DF_COMPLEXDF_H

#include <src/df/df.h>
#include <src/df/complexdf_base.h>
#include <src/df/complexdfinttask.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

class ComplexDFHalfDist;
class ComplexDFFullDist;

class ComplexDFDist : public DFDist, public ComplexDF_base {
  public:
    ComplexDFDist(const int nbas, const int naux, const std::array<std::shared_ptr<DFBlock>,2> block = std::array<std::shared_ptr<DFBlock>,2>{{nullptr, nullptr}},
                  std::shared_ptr<const ParallelDF> df = nullptr, std::shared_ptr<Matrix> data2 = nullptr);

    ComplexDFDist(const std::shared_ptr<const ParallelDF> df) : DFDist(df), ComplexDF_base() { }

    bool has_2index() const { return data2_.get() != nullptr; }
    size_t nbasis0() const { return nindex2_; }
    size_t nbasis1() const { return nindex1_; }
    size_t naux() const { return naux_; }

    // compute half transforms; c is dimensioned by nbasis_;
    std::shared_ptr<ComplexDFHalfDist> complex_compute_half_transform(const ZMatView c) const;

    // compute half transform using the third index. You get DFHalfDist with gamma/i/s (i.e., index are reordered)
    std::shared_ptr<ComplexDFHalfDist> complex_compute_half_transform_swap(const ZMatView c) const;

    // split up smalleri integrals into 6 dfdist objects
    std::vector<std::shared_ptr<const DFDist>> split_blocks() const override;

    // split up real and imaginary parts into 2 dfdist objects
    std::array<std::shared_ptr<const DFDist>,2> split_real_imag() const;

    // compute a J operator, given density matrices in AO basis
    std::shared_ptr<ZMatrix> complex_compute_Jop(const std::shared_ptr<const ZMatrix> den) const;
    std::shared_ptr<ZMatrix> complex_compute_Jop(const std::shared_ptr<const ComplexDF_base> o, const std::shared_ptr<const ZMatrix> den, const bool onlyonce = false) const;
    std::shared_ptr<ZMatrix> complex_compute_Jop_from_cd(std::shared_ptr<const ZVectorB> cd) const;
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
      for (int i = 0; i != 2*TBatch::Nblocks(); ++i) {
        blk[i] = block_[i];
      }

      int j2 = 0;
      for (auto& i2 : b2shell) {
        int j1 = 0;
        for (auto& i1 : b1shell) {
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
                       const double thr, const bool inverse, const double dum, const bool average = false, const std::shared_ptr<Matrix> data2 = nullptr) : ComplexDFDist(nbas, naux) {

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

      for (int i = 0; i != 2*TBatch::Nblocks(); ++i) {
        block_.push_back(std::make_shared<DFBlock>(adist_shell, adist_averaged, asize, b1size, b2size, astart, 0, 0));
      }

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

      assign_complex_blocks(*this);
    }

};

// TODO Functions not needed for RHF have not been implemented (transform_second, back_transform, etc.)
class ComplexDFHalfDist : public DFHalfDist, public ComplexDF_base {
  public:
    ComplexDFHalfDist(const std::shared_ptr<const ParallelDF> df, const int nocc) : DFHalfDist(df, nocc), ComplexDF_base() { }

    std::shared_ptr<ZMatrix> complex_form_2index(std::shared_ptr<const ComplexDFHalfDist> o, const double a, const bool swap = false) const;
    std::shared_ptr<ComplexDFHalfDist> complex_apply_J() const { return complex_apply_J(df_->data2()); }
    std::shared_ptr<ComplexDFHalfDist> complex_apply_J(const std::shared_ptr<const Matrix> o) const;

    std::shared_ptr<ComplexDFFullDist> complex_compute_second_transform(const ZMatView c) const;
    template<typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
    std::shared_ptr<ComplexDFFullDist> complex_compute_second_transform(std::shared_ptr<T> c) const { return compute_second_transform(*c); }

};


class ComplexDFFullDist : public DFFullDist, public ComplexDF_base {
  public:
    ComplexDFFullDist(const std::shared_ptr<const ParallelDF> df, const int nocc, const int nocc2) : DFFullDist(df, nocc, nocc2), ComplexDF_base() { }

    std::shared_ptr<ZMatrix> complex_form_4index(std::shared_ptr<const ComplexDFFullDist> o, const double a) const;
};

}

#endif
