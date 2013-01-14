//
// BAGEL - Parallel electron correlation program.
// Filename: df.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#ifndef __SRC_DF_DF_H
#define __SRC_DF_DF_H

#include <vector>
#include <memory>
#include <list>
#include <map>
#include <stdexcept>
#include <stddef.h>
#include <src/df/dfblock.h>
#include <src/scf/atom.h>
#include <src/util/matrix.h>
#include <src/df/dfinttask_old.h>
#include <src/df/dfinttask.h>
#include <src/scf/shell.h>
#include <src/util/timer.h>
#include <src/util/taskqueue.h>
#include <src/rysint/rysint.h>

namespace bagel {

class DFHalfDist;
class DFFullDist;

class ParallelDF : public std::enable_shared_from_this<ParallelDF> {
  protected:
    // blocks that this process has
    std::vector<std::shared_ptr<DFBlock> > block_;
    // hash key and process number
    std::map<int, int> global_table_;
    std::vector<std::pair<size_t, size_t> > atable_;

    // naux runs fastest, nindex2 runs slowest
    const size_t naux_;
    const size_t nindex1_;
    const size_t nindex2_;

    std::shared_ptr<const ParallelDF> df_;
    // data2_ is usually empty (except for the original DFDist)
    // AO two-index integrals ^ -1/2
    std::shared_ptr<Matrix> data2_;


  public:
    ParallelDF(const size_t, const size_t, const size_t);

    size_t naux() const { return naux_; }
    size_t nindex1() const { return nindex1_; }
    size_t nindex2() const { return nindex2_; }
    size_t size() const { return naux_*nindex1_*nindex2_; }

    std::vector<std::shared_ptr<DFBlock> >& block() { return block_; }
    const std::vector<std::shared_ptr<DFBlock> >& block() const { return block_; }
    std::shared_ptr<DFBlock> block(const size_t i) { return block_[i]; }
    std::shared_ptr<const DFBlock> block(const size_t i) const { return block_[i]; }

    void add_block(std::shared_ptr<DFBlock> o);

    std::shared_ptr<Matrix> form_2index(std::shared_ptr<const ParallelDF> o, const double a, const bool swap = false) const;
    std::unique_ptr<double[]> form_4index(std::shared_ptr<const ParallelDF> o, const double a, const bool swap = false) const;
    std::shared_ptr<Matrix> form_aux_2index(std::shared_ptr<const ParallelDF> o, const double a) const; 

    void daxpy(const double a, const std::shared_ptr<const ParallelDF> o);
    void scale(const double a);

    std::unique_ptr<double[]> get_block(const int i, const int id, const int j, const int jd, const int k, const int kd) const;

    const std::pair<size_t, size_t> atable(const int i) const { return atable_[i]; }
    const std::vector<std::pair<size_t, size_t> >& atable() const { return atable_; }

    const std::shared_ptr<const ParallelDF> df() const { return df_; }
    std::shared_ptr<const Matrix> data2() const { return data2_; }

    int get_node(const int shelloffset) const;

    // compute a J operator, given density matrices in AO basis
    std::shared_ptr<Matrix> compute_Jop(const std::shared_ptr<const Matrix> den) const;
    std::shared_ptr<Matrix> compute_Jop(const std::shared_ptr<const ParallelDF> o, const std::shared_ptr<const Matrix> den) const;

    std::unique_ptr<double[]> compute_cd(const std::shared_ptr<const Matrix> den, std::shared_ptr<const Matrix> dat2 = std::shared_ptr<const Matrix>()) const;

};


class DFDist : public ParallelDF {
  friend class DFIntTask_OLD<DFDist>;
  friend class DFFullDist;
  friend class DFHalfDist;
  protected:
    std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input);

    virtual void common_init1(const std::vector<std::shared_ptr<const Atom> >&,
                              const std::vector<std::shared_ptr<const Atom> >&,
                              const std::vector<std::shared_ptr<const Atom> >&, const double thresh, const bool compute_inv) { assert(false); }
    void common_init2(const std::vector<std::shared_ptr<const Atom> >&,
                      const std::vector<std::shared_ptr<const Atom> >&,
                      const std::vector<std::shared_ptr<const Atom> >&, const double thresh, const bool compute_inv);
    void make_table(const int nmax);
    std::tuple<int, std::vector<std::shared_ptr<const Shell> > > get_ashell(const std::vector<std::shared_ptr<const Shell> >& all) const;

  public:
    // construction of a block from AO integrals
    DFDist(const int nbas, const int naux, const std::vector<std::shared_ptr<const Atom> >& atoms,
                                           const std::vector<std::shared_ptr<const Atom> >& aux_atoms, const double thr, const bool inverse, const double dum)
      : ParallelDF(naux, nbas, nbas) {
    }

    DFDist(const std::shared_ptr<const ParallelDF> df) : ParallelDF(df->naux(), df->nindex1(), df->nindex2()) { df_ = df; }

    bool has_2index() const { return data2_.get() != nullptr; };
    size_t nbasis0() const { return nindex2_; };
    size_t nbasis1() const { return nindex1_; };
    size_t naux() const { return naux_; };

    void add_direct_product(std::vector<const double*> a, std::vector<const double*> b, const double fac);

    // compute half transforms; c is dimensioned by nbasis_;
    std::shared_ptr<DFHalfDist> compute_half_transform(const double* c, const size_t nocc) const;
    std::shared_ptr<DFHalfDist> compute_half_transform(const std::shared_ptr<const Matrix> c) const { return compute_half_transform(c->data(), c->mdim()); }

    // compute half transform using the third index. You get DFHalfDist with gamma/i/s (i.e., index are reordered)
    std::shared_ptr<DFHalfDist> compute_half_transform_swap(const double* c, const size_t nocc) const;
    std::shared_ptr<DFHalfDist> compute_half_transform_swap(const std::shared_ptr<const Matrix> c) const { return compute_half_transform_swap(c->data(), c->mdim()); }

};


template<class TBatch>
class DFDist_ints : public DFDist {
  protected:
    void common_init1(const std::vector<std::shared_ptr<const Atom> >& atoms0,
                      const std::vector<std::shared_ptr<const Atom> >& atoms1,
                      const std::vector<std::shared_ptr<const Atom> >& aux_atoms, const double thresh, const bool compute_inv) {
      Timer time;

      // 3index Integral is now made in DFBlock.
      std::vector<std::shared_ptr<const Shell> > ashell, b1shell, b2shell;
      for (auto& i : aux_atoms) ashell.insert(ashell.end(), i->shells().begin(), i->shells().end());
      for (auto& i : atoms1) b1shell.insert(b1shell.end(), i->shells().begin(), i->shells().end());
      for (auto& i : atoms0) b2shell.insert(b2shell.end(), i->shells().begin(), i->shells().end());

      // get number of dfblocks
      const int nblocks = TBatch::nblocks();
      block_.resize(nblocks);

      int astart;
      std::vector<std::shared_ptr<const Shell> > myashell;
      std::tie(astart, myashell) = get_ashell(ashell);

      // make empty dfblocks
      for (int i = 0; i != nblocks; ++i)
        block_[i] = std::shared_ptr<DFBlock>(new DFBlock(myashell, b1shell, b2shell, astart, 0, 0));

      // making a task list
      std::vector<DFIntTask<TBatch> > tasks;
      tasks.reserve(b1shell.size()*b2shell.size()*myashell.size());

      // TODO this is not general, but for the time being I plan to have full basis functions (and limited aux basis functions); 
      assert(b1shell == b2shell);

      const std::shared_ptr<const Shell> i3(new Shell(myashell.front()->spherical()));

      auto j2 = block_[0]->b2off().begin();
      for (auto& i2 : b2shell) {
        auto j1 = block_[0]->b1off().begin();
        for (auto& i1 : b1shell) {
          // TODO using symmetry. This assumes that swap(i1, i2) integrals are also located in this block, which might not be the case in general.
          if (*j1 <= *j2) {
            auto j0 = block_[0]->aoff().begin();
            for (auto& i0 : myashell) {
              tasks.push_back(DFIntTask<TBatch>(std::array<std::shared_ptr<const Shell>,4>{{i3, i0, i1, i2}}, std::vector<int>{*j2, *j1, *j0}, block_));
              ++j0;
            }
          }
          ++j1;
        }
        ++j2;
      }
      TaskQueue<DFIntTask<TBatch> > tq(tasks);
      tq.compute(resources__->max_num_threads());
      time.tick_print("3-index ints");
    }

  public:
    DFDist_ints(const int nbas, const int naux, const std::vector<std::shared_ptr<const Atom> >& atoms,
                                           const std::vector<std::shared_ptr<const Atom> >& aux_atoms, const double thr, const bool inverse, const double dum)
      : DFDist(nbas, naux, atoms, aux_atoms, thr, inverse, dum) {
      common_init1(atoms, atoms, aux_atoms, thr, inverse);
      common_init2(atoms, atoms, aux_atoms, thr, inverse);
    }

};


class DFHalfDist : public ParallelDF {
  friend class DFFullDist;
  protected:
    std::shared_ptr<DFHalfDist> apply_J(const std::shared_ptr<const Matrix> o) const;

  public:
    DFHalfDist(const std::shared_ptr<const ParallelDF> df, const int nocc) : ParallelDF(df->naux(), nocc, df->nindex2()) { df_ = df; }

    size_t nocc() const { return nindex1_; };
    size_t nbasis() const { return nindex2_; };

    std::shared_ptr<DFFullDist> compute_second_transform(const double* c, const size_t nocc) const;
    std::shared_ptr<DFFullDist> compute_second_transform(const std::shared_ptr<const Matrix> c) const { return compute_second_transform(c->data(), c->mdim()); }
    std::shared_ptr<DFDist> back_transform(const double* c) const;
    std::shared_ptr<DFDist> back_transform(const std::shared_ptr<const Matrix> c) const { assert(c->mdim() == nindex1_); return back_transform(c->data()); }

    std::shared_ptr<DFHalfDist> copy() const; 
    std::shared_ptr<DFHalfDist> clone() const; 

    void rotate_occ(const double* d);
    std::shared_ptr<DFHalfDist> apply_density(const double* d) const;

    std::shared_ptr<Matrix> compute_Kop_1occ(const double* den, const double a) const;

    std::shared_ptr<DFHalfDist> apply_J() const { return apply_J(df_->data2()); }
    std::shared_ptr<DFHalfDist> apply_JJ() const { return apply_J(std::shared_ptr<Matrix>(new Matrix(*df_->data2()**df_->data2()))); }
    std::shared_ptr<DFHalfDist> apply_J(const std::shared_ptr<const DFDist> d) const { return apply_J(d->data2()); }
    std::shared_ptr<DFHalfDist> apply_JJ(const std::shared_ptr<const DFDist> d) const { return apply_J(std::shared_ptr<Matrix>(new Matrix(*d->data2()**d->data2()))); }

};


class DFFullDist : public ParallelDF {
  protected:
    std::shared_ptr<DFFullDist> apply_J(const std::shared_ptr<const Matrix> o) const;

  public:
    DFFullDist(const std::shared_ptr<const ParallelDF> df, const int nocc1, const int nocc2) : ParallelDF(df->naux(), nocc1, nocc2) { df_ = df; }

    int nocc1() const { return nindex1_; }
    int nocc2() const { return nindex2_; }

    std::shared_ptr<DFFullDist> copy() const; 
    std::shared_ptr<DFFullDist> clone() const; 

    std::shared_ptr<DFHalfDist> back_transform(const double* c) const;
    std::shared_ptr<DFHalfDist> back_transform(const std::shared_ptr<const Matrix> c) const { assert(c->mdim() == nindex2_); return back_transform(c->data()); }

    void symmetrize();

    // 2RDM contractions
    // special function for RHF
    std::shared_ptr<DFFullDist> apply_closed_2RDM() const;
    // special function for UHF
    std::shared_ptr<DFFullDist> apply_uhf_2RDM(const double* rdma, const double* rdmb) const;
    // general case with closed orbitals
    std::shared_ptr<DFFullDist> apply_2rdm(const double* rdm, const double* rdm1, const int nclosed, const int nact) const;
    // general case without closed orbitals
    std::shared_ptr<DFFullDist> apply_2rdm(const double* rdm) const;

    std::unique_ptr<double[]> form_4index_1fixed(const std::shared_ptr<const DFFullDist> o, const double a, const size_t n) const;

    // utility functions
    std::shared_ptr<Matrix> form_aux_2index_apply_J(const std::shared_ptr<const DFFullDist> o, const double a) const;
    void set_product(const std::shared_ptr<const DFFullDist>, const std::unique_ptr<double[]>&, const int jdim, const size_t offset); 

    std::shared_ptr<DFFullDist> apply_J() const { return apply_J(df_->data2()); }
    std::shared_ptr<DFFullDist> apply_JJ() const { return apply_J(std::shared_ptr<Matrix>(new Matrix(*df_->data2()**df_->data2()))); }
    std::shared_ptr<DFFullDist> apply_J(const std::shared_ptr<const ParallelDF> d) const { return apply_J(d->data2()); }
    std::shared_ptr<DFFullDist> apply_JJ(const std::shared_ptr<const ParallelDF> d) const { return apply_J(std::shared_ptr<Matrix>(new Matrix(*d->data2()**d->data2()))); }

};

}


#endif
