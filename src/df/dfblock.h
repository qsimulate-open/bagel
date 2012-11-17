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
    size_t asize_;
    size_t b1size_;
    size_t b2size_;

    // a set of offsets of this block in the entire DF integrals
    size_t astart_; 
    size_t b1start_; 
    size_t b2start_; 

    void common_init();

    std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input);

  public:
    DFBlock(std::vector<std::shared_ptr<const Shell> > a,
            std::vector<std::shared_ptr<const Shell> > b1,
            std::vector<std::shared_ptr<const Shell> > b2,
            const int as, const int b1s, const int b2s)
     : aux_(a), b1_(b1), b2_(b2), asize_(0lu), b1size_(0lu), b2size_(0lu), astart_(as), b1start_(b1s), b2start_(b2s) {

      for (auto& i : aux_) { aoff_.push_back(asize_); asize_ += i->nbasis(); }
      for (auto& i : b1_)  { b1off_.push_back(b1size_); b1size_ += i->nbasis(); }
      for (auto& i : b2_)  { b2off_.push_back(b2size_); b2size_ += i->nbasis(); }

      common_init();
    };


    DFBlock(std::unique_ptr<double[]>& d) : data_(std::move(d)) { };

    // direct access to data will be disabled once implementation is done
    double* get() { return data_.get(); };
    const double* get() const { return data_.get(); };
};

}


#endif
