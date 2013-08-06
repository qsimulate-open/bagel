//
// BAGEL - Parallel electron correlation program.
// Filename: reldfhalf.h
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


#ifndef __SRC_REL_RELDFHALF_H
#define __SRC_REL_RELDFHALF_H

#include <memory>
#include <string>
#include <map>
#include <src/wfn/reference.h>
#include <src/math/zmatrix.h>
#include <src/math/matrix.h>
#include <src/rel/alpha.h>
#include <src/rel/breit2index.h>
#include <src/rel/reldf.h>
#include <src/df/df.h>

namespace bagel {

class RelDF;

class RelDFHalf : public RelDFBase {
  protected:
    std::array<std::shared_ptr<DFHalfDist>,2> dfhalf_;
    std::array<std::shared_ptr<DFHalfDist>,2> df2_;
    bool split_;

  public:
    RelDFHalf(std::shared_ptr<const RelDF>, std::vector<std::shared_ptr<const SpinorInfo>> bas,
                  std::array<std::shared_ptr<const Matrix>,4>, std::array<std::shared_ptr<const Matrix>,4>);

    RelDFHalf(std::array<std::shared_ptr<DFHalfDist>,2> data, std::pair<int,int> cartesian, std::vector<std::shared_ptr<const SpinorInfo>> bas);
    RelDFHalf(const RelDFHalf& o);

    std::array<std::shared_ptr<DFHalfDist>, 2> get_data() const { return dfhalf_; }
    std::shared_ptr<DFHalfDist> get_real() const { return dfhalf_[0]; }
    std::shared_ptr<DFHalfDist> get_imag() const { return dfhalf_[1]; }

    bool matches(std::shared_ptr<const RelDFHalf>) const;
    bool alpha_matches(std::shared_ptr<const Breit2Index>) const;
    bool alpha_matches(std::shared_ptr<const RelDFHalf>) const;
    std::shared_ptr<RelDFHalf> multiply_breit2index(std::shared_ptr<const Breit2Index> b2i) const;

    std::shared_ptr<RelDFHalf> copy() const { return std::make_shared<RelDFHalf>(*this); }
    std::shared_ptr<RelDFHalf> apply_J() const;

    // zaxpy
    void zaxpy(std::complex<double> a, std::shared_ptr<const RelDFHalf> o);

    // for the zgemm3m-like algorithm
    void set_sum_diff();
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
