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
#include <src/scf/shell.h>
#include <src/rysint/rysint.h>

namespace bagel {

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

    void ao_init();

    std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input);

  public:
    // construction of a block from AO integrals
    DFBlock(std::vector<std::shared_ptr<const Shell> > a, std::vector<std::shared_ptr<const Shell> > b1, std::vector<std::shared_ptr<const Shell> > b2,
            const int as, const int b1s, const int b2s);

    // construction of a block from data (unique_ptr<double[]>)
    DFBlock(std::unique_ptr<double[]>& d, const size_t a, const size_t b1, const size_t b2,
                                          const int as, const int b1s, const int b2s)
     : data_(std::move(d)), asize_(a), b1size_(b1), b2size_(b2), astart_(as), b1start_(b1s), b2start_(b2s) { };

    std::shared_ptr<DFBlock> transform_second(const double* const c, const int nocc) const;
    std::shared_ptr<DFBlock> transform_third(const double* const c, const int nocc) const;

    std::shared_ptr<DFBlock> clone() const;
    std::shared_ptr<DFBlock> copy() const;

    // dimensions of the block
    size_t asize() const { return asize_; }
    size_t b1size() const { return b1size_; }
    size_t b2size() const { return b2size_; }

    size_t size() const { return asize_*b1size_*b2size_; };

    // a set of offsets of this block in the entire DF integrals
    size_t astart() const { return astart_; }
    size_t b1start() const { return b1start_; }
    size_t b2start() const { return b2start_; }

    // some math functions
    DFBlock& operator+=(const DFBlock& o);
    DFBlock& operator-=(const DFBlock& o);
    void daxpy(const double a, const DFBlock& o);
    void daxpy(const double a, const std::shared_ptr<const DFBlock> o) { daxpy(a, *o); }
    void scale(const double a);

    // some additional functions
    // symmetrize b1 and b2 (assuming b1size_ == b2size_)
    void symmetrize();

    // TODO direct access to data will be disabled once implementation is done
    double* get() { return data_.get(); }
    const double* get() const { return data_.get(); }
    double& operator[](const size_t i) { return data_[i]; }
};

}


#endif
