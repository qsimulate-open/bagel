//
// BAGEL - Parallel electron correlation program.
// Filename: dfblockt.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_DF_DFBLOCKT_H
#define __SRC_DF_DFBLOCKT_H

#include <src/df/df.h>

namespace bagel {

/*
    DFDistT is a slice of 3-index DF integrals. Distributed by the second and third index.
    Note that the date is tranposed!
*/

class DFDistT {
  protected:
    std::vector<std::shared_ptr<Matrix>> data_;

    // first dimension is naux_ (global)
    const size_t naux_;
    // original dimensions
    const size_t nindex1_;
    const size_t nindex2_;

    // distribution information
    std::shared_ptr<const StaticDist> dist_;

    // second and third dimension
    size_t start_;
    size_t size_;

    std::shared_ptr<const ParallelDF> df_;

  public:
    // CAUTION this constructor should be called **COLLECTIVELY**!! Otherwise the program hangs.
    DFDistT(std::shared_ptr<const ParallelDF> in, std::shared_ptr<const StaticDist> dist = std::shared_ptr<const StaticDist>());

    DFDistT(const size_t naux, std::shared_ptr<const StaticDist> dist, const size_t n1, const size_t n2,
            const std::shared_ptr<const ParallelDF>);

    std::shared_ptr<DFDistT> clone() const;
    std::shared_ptr<DFDistT> apply_J(std::shared_ptr<const Matrix> d) const;
    std::shared_ptr<DFDistT> apply_J() const { return apply_J(df_->data2()); }
    std::vector<std::shared_ptr<Matrix>> form_aux_2index(std::shared_ptr<const DFDistT> o, const double a) const;

    void get_paralleldf(std::shared_ptr<ParallelDF>) const;

    int size() const { return size_; }
    int start() const { return start_; }
    int nblocks() const { return data_.size(); }
    const double* data() const { assert(data_.size() == 1); return data(0); }
    const double* data(const int i) const { return data_[i]->data(); }

    std::vector<std::shared_ptr<Matrix>> get_slice(const int start, const int end) const;

    std::shared_ptr<const ParallelDF> df() const { return df_; }
    void discard_df() { df_.reset(); }
};

}

#endif
