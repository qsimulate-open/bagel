//
// BAGEL - Parallel electron correlation program.
// Filename: reldfhalf_london.h
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


#ifndef __SRC_LONDON_RELDFHALF_LONDON_H
#define __SRC_LONDON_RELDFHALF_LONDON_H

#include <memory>
#include <string>
#include <map>
#include <src/wfn/reference.h>
#include <src/math/zmatrix.h>
#include <src/math/matrix.h>
#include <src/rel/alpha.h>
#include <src/rel/breit2index.h>
#include <src/london/reldf_london.h>
#include <src/df/df_london.h>

namespace bagel {

class RelDF_London;

class RelDFHalf_London : public RelDFBase {
  protected:
    std::array<std::shared_ptr<DFHalfDist_London>,2> dfhalf_;
    std::array<std::shared_ptr<DFHalfDist_London>,2> df2_;
    bool split_;

  public:
    RelDFHalf_London(std::shared_ptr<const RelDF_London>, std::vector<std::shared_ptr<const SpinorInfo>> bas,
                     std::array<std::shared_ptr<const ZMatrix>,4>, std::array<std::shared_ptr<const ZMatrix>,4>);

    RelDFHalf_London(std::array<std::shared_ptr<DFHalfDist_London>,2> data, std::pair<int,int> cartesian, std::vector<std::shared_ptr<const SpinorInfo>> bas);
    RelDFHalf_London(const RelDFHalf_London& o);

    std::array<std::shared_ptr<DFHalfDist_London>, 2> get_data() const { return dfhalf_; }
    std::shared_ptr<DFHalfDist_London> get_real() const { return dfhalf_[0]; }
    std::shared_ptr<DFHalfDist_London> get_imag() const { return dfhalf_[1]; }

    bool matches(std::shared_ptr<const RelDFHalf_London>) const;
    bool alpha_matches(std::shared_ptr<const Breit2Index>) const;
    bool alpha_matches(std::shared_ptr<const RelDFHalf_London>) const;
    std::shared_ptr<RelDFHalf_London> multiply_breit2index(std::shared_ptr<const Breit2Index> b2i) const;

    std::shared_ptr<RelDFHalf_London> copy() const { return std::make_shared<RelDFHalf_London>(*this); }
    std::shared_ptr<RelDFHalf_London> apply_J() const;

    // zaxpy
    void ax_plus_y(std::complex<double> a, std::shared_ptr<const RelDFHalf_London> o);

    // for the zgemm3m-like algorithm
    void set_sum_diff();
    void discard_sum_diff() { df2_ = std::array<std::shared_ptr<DFHalfDist_London>,2>(); }
    std::shared_ptr<DFHalfDist_London> sum() const { return df2_[0]; }
    std::shared_ptr<DFHalfDist_London> diff() const { return df2_[1]; }

    std::complex<double> fac() const { assert(basis_.size() == 1); return basis_[0]->fac(cartesian_); }
    std::list<std::shared_ptr<RelDFHalf_London>> split(const bool docopy = false);
    bool split_status() const { return split_; }

};


// Half-transformed DF objects when backtransforming to AO
class RelDFHalfB_London {
  protected:
    std::array<std::shared_ptr<DFHalfDist_London>,2> dfhalf_;
    const int basis_;
    const int alpha_;
  public:
    RelDFHalfB_London(std::array<std::shared_ptr<DFHalfDist_London>,2> data, const int basis, const int alpha) : dfhalf_(data), basis_(basis), alpha_(alpha) { }

    int basis() const { return basis_; }
    std::shared_ptr<DFDist_London> back_transform(std::shared_ptr<const ZMatrix>, std::shared_ptr<const ZMatrix>, const bool imag = false) const;
};

}

#endif
