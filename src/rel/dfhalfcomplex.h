//
// BAGEL - Parallel electron correlation program.
// Filename: dfhalfcomplex.h
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


#ifndef __SRC_REL_DFHALFCOMPLEX_H
#define __SRC_REL_DFHALFCOMPLEX_H

#include <memory>
#include <string>
#include <map>
#include <src/wfn/reference.h>
#include <src/util/zmatrix.h>
#include <src/util/matrix.h>
#include <src/rel/dfdata.h>
#include <src/df/df.h>

namespace bagel {


class DFHalfComplex {
  protected:
    std::array<std::shared_ptr<DFHalfDist>, 2> dfhalf_;
    std::pair<const int, const int> coord_;
    std::pair<const int, const int> basis_; 
    int dim_;

    std::array<std::shared_ptr<DFHalfDist>, 2> df2_;

  public:
    DFHalfComplex(const std::shared_ptr<const DFDist>, std::shared_ptr<const Matrix>, std::shared_ptr<const Matrix>, 
                  const bool, std::pair<const int, const int>, std::pair<const int, const int>);

    std::array<std::shared_ptr<DFHalfDist>, 2> get_data() { return dfhalf_; }
    std::shared_ptr<DFHalfDist> get_real() { return dfhalf_[0]; }
    std::shared_ptr<DFHalfDist> get_imag() { return dfhalf_[1]; }
    std::pair<const int, const int> coord() { return coord_; }
    std::pair<const int, const int> basis() { return basis_; }

    const std::tuple<int, int, int, int> compute_index_Jop(std::pair<const int, const int>, std::pair<const int, const int>);
    const std::tuple<int, int, int, int> compute_index_Exop(std::pair<const int, const int>, std::pair<const int, const int>); 
    std::complex<double> compute_coeff(std::pair<const int, const int>, std::pair<const int, const int>);
    const int coeff_matrix() const;

    // zaxpy
    void zaxpy(std::complex<double> a, std::shared_ptr<const DFHalfComplex> o);

    // for the zgemm3m-like algorithm
    void set_sum_diff();
    std::shared_ptr<DFHalfDist> sum() const { return df2_[0]; } 
    std::shared_ptr<DFHalfDist> diff() const { return df2_[1]; } 
};

}

#endif
