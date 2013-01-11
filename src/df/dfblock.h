//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock.h
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

#ifndef __SRC_DF_DFBLOCK_H
#define __SRC_DF_DFBLOCK_H

#include <numeric>
#include <memory>
#include <list>
#include <cassert>
#include <vector>
#include <src/util/matrix.h>
#include <src/scf/shell.h>
#include <src/util/timer.h>
#include <src/util/taskqueue.h>
#include <src/rysint/rysint.h>

namespace bagel {

/*
    DFBlock is a slice of 3-index DF integrals. Distributed by the first index
*/

class DFBlock {
  friend class DFIntTask;
  protected:
    // aux_ runs fastest, b2_ runs slowest
    std::unique_ptr<double[]> data_;
 
    std::vector<std::shared_ptr<const Shell> > aux_;
    std::vector<std::shared_ptr<const Shell> > b1_;
    std::vector<std::shared_ptr<const Shell> > b2_;
    // offset within the block
    std::vector<int> aoff_;
    std::vector<int> b1off_;
    std::vector<int> b2off_;

    // dimensions of the block
    const size_t asize_;
    const size_t b1size_;
    const size_t b2size_;

    // a set of offsets of this block in the entire DF integrals
    const size_t astart_; 
    const size_t b1start_; 
    const size_t b2start_; 

    virtual std::shared_ptr<Integral> compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) {
      assert(false);
      return std::shared_ptr<Integral>();
    }

  public:
    // construction of a block from AO integrals
    DFBlock(std::vector<std::shared_ptr<const Shell> > a, std::vector<std::shared_ptr<const Shell> > b1, std::vector<std::shared_ptr<const Shell> > b2,
            const int as, const int b1s, const int b2s);

    // construction of a block from data (unique_ptr<double[]>)
    DFBlock(std::unique_ptr<double[]>& d, const size_t a, const size_t b1, const size_t b2,
                                          const int as, const int b1s, const int b2s)
     : data_(std::move(d)), asize_(a), b1size_(b1), b2size_(b2), astart_(as), b1start_(b1s), b2start_(b2s) { }

    std::shared_ptr<DFBlock> transform_second(const double* const c, const int nocc, const bool trans = false) const;
    std::shared_ptr<DFBlock> transform_third(const double* const c, const int nocc, const bool trans = false) const;

    std::shared_ptr<DFBlock> clone() const;
    std::shared_ptr<DFBlock> copy() const;
    void zero() { std::fill_n(data_.get(), size(), 0.0); }

    // dimensions of the block
    size_t asize() const { return asize_; }
    size_t b1size() const { return b1size_; }
    size_t b2size() const { return b2size_; }

    size_t size() const { return asize_*b1size_*b2size_; }

    // a set of offsets of this block in the entire DF integrals
    size_t astart() const { return astart_; }
    size_t b1start() const { return b1start_; }
    size_t b2start() const { return b2start_; }

    // TODO direct access to data will be disabled once implementation is done
    double* get() { return data_.get(); }
    const double* get() const { return data_.get(); }
    double& operator[](const size_t i) { return data_[i]; }

    // some math functions
    DFBlock& operator+=(const DFBlock& o);
    DFBlock& operator-=(const DFBlock& o);
    void daxpy(const double a, const DFBlock& o);
    void daxpy(const double a, const std::shared_ptr<const DFBlock> o) { daxpy(a, *o); }
    void scale(const double a);

    // add ab^+  to this.
    void add_direct_product(const double* a, const double* b, const double fac);

    // some additional functions
    // symmetrize b1 and b2 (assuming b1size_ == b2size_)
    void symmetrize();

    // 2RDM contractions
    std::shared_ptr<DFBlock> apply_rhf_2RDM() const; 
    std::shared_ptr<DFBlock> apply_uhf_2RDM(const double*, const double*) const; 
    std::shared_ptr<DFBlock> apply_2RDM(const double* rdm, const double* rdm1, const int nclosed, const int nact) const;
    std::shared_ptr<DFBlock> apply_2RDM(const double* rdm) const;

    // Form 2- and 4-index integrals
    std::shared_ptr<Matrix> form_2index(const std::shared_ptr<const DFBlock> o, const double a) const;
    std::unique_ptr<double[]> form_4index(const std::shared_ptr<const DFBlock> o, const double a) const;
    // slowest index of o is fixed to n
    std::unique_ptr<double[]> form_4index_1fixed(const std::shared_ptr<const DFBlock> o, const double a, const size_t n) const;
    std::shared_ptr<Matrix> form_aux_2index(const std::shared_ptr<const DFBlock> o, const double a) const;

    std::unique_ptr<double[]> form_vec(const std::shared_ptr<const Matrix> den) const;
    std::shared_ptr<Matrix> form_mat(const double* fit) const;

    void contrib_apply_J(const std::shared_ptr<const DFBlock> o, const std::shared_ptr<const Matrix> mat);

    void copy_block(const std::unique_ptr<double[]>& o, const int jdim, const size_t offset);
    // compute (D|ia)(ia|j) and set to the location specified by the offset
    std::unique_ptr<double[]> form_Dj(const std::unique_ptr<double[]>& o, const int jdim) const;

    // CAUTION, ist, jst, and kst are absolute number (NOT relative to astart_, ...). Returns double[] whose size is i*j*k 
    std::unique_ptr<double[]> get_block(const int ist, const int i, const int jst, const int j, const int kst, const int k) const; 

    // use with caution
    std::unique_ptr<double[]> release_data() { return std::move(data_); }
};

}

#include <src/df/dfinttask.h>

namespace bagel {

template<class TBatch>
class DFBlock_ints : public DFBlock {
  protected:
    std::shared_ptr<Integral> compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) {
      std::shared_ptr<TBatch> eribatch(new TBatch(input, 2.0));
      eribatch->compute();
      return eribatch;
    }

    void ao_init() {
      // allocation of the data area
      data_ = std::unique_ptr<double[]>(new double[asize_*b1size_*b2size_]);

      const std::shared_ptr<const Shell> i3(new Shell(aux_.front()->spherical()));

      // making a task list
      std::vector<DFIntTask> tasks;
      tasks.reserve(b1_.size()*b2_.size()*aux_.size());

      // TODO this is not general, but for the time being I plan to have full basis functions (and limited aux basis functions); 
      assert(b1_ == b2_);
  
      auto j2 = b2off_.begin();
      for (auto& i2 : b2_) { 
        auto j1 = b1off_.begin();
        for (auto& i1 : b1_) { 
          // TODO using symmetry. This assumes that swap(i1, i2) integrals are also located in this block, which might not be the case in general.
          if (*j1 <= *j2) {
            auto j0 = aoff_.begin();
            for (auto& i0 : aux_)
              tasks.push_back(DFIntTask(std::array<std::shared_ptr<const Shell>,4>{{i3, i0, i1, i2}}, std::vector<int>{*j2, *j1, *j0++}, this));
          }
          ++j1;
        }
        ++j2;
      }
      TaskQueue<DFIntTask> tq(tasks);
      tq.compute(resources__->max_num_threads());
    }

  public:
    // construction of a block from AO integrals
    DFBlock_ints(std::vector<std::shared_ptr<const Shell> > a, std::vector<std::shared_ptr<const Shell> > b1, std::vector<std::shared_ptr<const Shell> > b2,
             const int as, const int b1s, const int b2s) : DFBlock(a, b1, b2, as, b1s, b2s) { ao_init(); }
    
}; 


}


#endif
