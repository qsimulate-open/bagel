//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reldfhalf.h
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


#ifndef __SRC_DF_RELDFHALF_H
#define __SRC_DF_RELDFHALF_H

#include <map>
#include <src/wfn/reference.h>
#include <src/df/breit2index.h>
#include <src/df/reldf.h>

namespace bagel {

class RelDF;

class RelDFHalf : public RelDFBase {
  protected:
    std::array<std::shared_ptr<DFHalfDist>,2> dfhalf_;
    bool split_;

    // for efficiency... blame use of mutable
    mutable std::array<std::shared_ptr<DFHalfDist>,2> df2_;

  public:
    RelDFHalf(std::shared_ptr<const RelDF>, std::vector<std::shared_ptr<const SpinorInfo>> bas,
                  std::array<std::shared_ptr<const Matrix>,4>, std::array<std::shared_ptr<const Matrix>,4>);

    RelDFHalf(std::array<std::shared_ptr<DFHalfDist>,2> data, std::pair<int,int> cartesian, std::vector<std::shared_ptr<const SpinorInfo>> bas);
    RelDFHalf(const RelDFHalf& o);

    int nocc() const { return dfhalf_[0]->nocc(); }

    std::array<std::shared_ptr<DFHalfDist>, 2> get_data() const { return dfhalf_; }
    std::shared_ptr<DFHalfDist> get_real() const { return dfhalf_[0]; }
    std::shared_ptr<DFHalfDist> get_imag() const { return dfhalf_[1]; }

    bool matches(std::shared_ptr<const RelDFHalf>) const;
    bool alpha_matches(std::shared_ptr<const Breit2Index>) const;
    bool alpha_matches(std::shared_ptr<const RelDFHalf>) const;
    std::shared_ptr<RelDFHalf> multiply_breit2index(std::shared_ptr<const Breit2Index> b2i) const;

    std::shared_ptr<RelDFHalf> copy() const { return std::make_shared<RelDFHalf>(*this); }
    std::shared_ptr<RelDFHalf> apply_J() const;
    std::shared_ptr<RelDFHalf> apply_JJ() const;

    std::shared_ptr<RelDFHalf> merge_b1(std::shared_ptr<RelDFHalf> o) const;
    std::shared_ptr<RelDFHalf> slice_b1(const int slice_start, const int slice_size) const;

    void ax_plus_y(std::complex<double> a, std::shared_ptr<const RelDFHalf> o);
    std::shared_ptr<RelDFHalf> transform_occ(std::shared_ptr<const ZMatrix> rdm1) const;

    // for the zgemm3m-like algorithm
    void set_sum_diff() const;
    void discard_sum_diff() const { df2_ = std::array<std::shared_ptr<DFHalfDist>,2>(); }
    std::shared_ptr<DFHalfDist> sum() const { return df2_[0]; }
    std::shared_ptr<DFHalfDist> diff() const { return df2_[1]; }

    std::complex<double> fac() const { assert(basis_.size() == 1); return basis_[0]->fac(cartesian_); }
    std::list<std::shared_ptr<RelDFHalf>> split(const bool docopy = false);
    bool split_status() const { return split_; }

};


// Half-transformed DF objects when backtransforming to AO
class RelDFHalfB {
  protected:
    std::array<std::shared_ptr<DFHalfDist>,2> dfhalf_;
    const int basis_;
    const int alpha_;
  public:
    RelDFHalfB(std::array<std::shared_ptr<DFHalfDist>,2> data, const int basis, const int alpha) : dfhalf_(data), basis_(basis), alpha_(alpha) { }

    int basis() const { return basis_; }
    std::shared_ptr<DFDist> back_transform(std::shared_ptr<const Matrix>, std::shared_ptr<const Matrix>, const bool imag = false) const;
};

}

#endif
