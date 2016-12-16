//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dfdistt.h
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


#ifndef __SRC_DF_DFDISTT_H
#define __SRC_DF_DFDISTT_H

#include <src/df/df.h>

namespace bagel {

/*
    DFDistT is a 3-index DF integral object, which is distributed by the second and third indices.
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
    size_t bstart_;
    size_t bsize_;

    std::shared_ptr<const ParallelDF> df_;

  public:
    // CAUTION this constructor should be called **COLLECTIVELY**!! Otherwise the program hangs.
    DFDistT(std::shared_ptr<const ParallelDF> in, std::shared_ptr<const StaticDist> dist = nullptr);

    DFDistT(const size_t naux, std::shared_ptr<const StaticDist> dist, const size_t n1, const size_t n2,
            const std::shared_ptr<const ParallelDF>);

    std::shared_ptr<DFDistT> clone() const;
    std::shared_ptr<DFDistT> apply_J(std::shared_ptr<const Matrix> d) const;
    std::shared_ptr<DFDistT> apply_J() const { return apply_J(df_->data2()); }
    std::vector<std::shared_ptr<Matrix>> form_aux_2index(std::shared_ptr<const DFDistT> o, const double a) const;

    void get_paralleldf(std::shared_ptr<ParallelDF>) const;

    size_t naux() const { return naux_; }
    size_t nindex1() const { return nindex1_; }
    size_t nindex2() const { return nindex2_; }

    int bsize() const { return bsize_; }
    int bstart() const { return bstart_; }
    int nblocks() const { return data_.size(); }
    double* data() { assert(data_.size() == 1); return data(0); }
    double* data(const int i) { return data_[i]->data(); }
    const double* data() const { assert(data_.size() == 1); return data(0); }
    const double* data(const int i) const { return data_[i]->data(); }

    // returns the process that has the data
    int locate(const size_t, const size_t n) const { return std::get<0>(dist_->locate(n)); }
    size_t offset(const size_t, const size_t n) const { return naux_*std::get<1>(dist_->locate(n)); }

    std::vector<std::shared_ptr<Matrix>> get_slice(const int start, const int end) const;

    std::shared_ptr<const ParallelDF> df() const { return df_; }
    void discard_df() { df_.reset(); }

    std::shared_ptr<Matrix> replicate(const int i = 0) const;
};

}

#endif
