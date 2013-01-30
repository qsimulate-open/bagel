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
#include <src/df/df.h>

namespace bagel {

class DFHalfComplex {
  protected:
    std::array<std::shared_ptr<DFHalfDist>, 2> dfdata_;
    std::array<std::shared_ptr<Matrix>, 2> data_;
    std::pair<const int, const int> coord_;
    std::pair<const int, const int> basis_; 
    int dim_;
    void initialize_data_();

  public:
    DFHalfComplex(const std::shared_ptr<const DFDist>, std::shared_ptr<const Matrix>, std::shared_ptr<const Matrix>, 
                  const bool, std::pair<const int, const int>, std::pair<const int, const int>);

    std::array<std::shared_ptr<DFHalfDist>, 2> get_data() { return dfdata_; }
    std::shared_ptr<DFHalfDist> get_real() { return dfdata_[0]; }
    std::shared_ptr<DFHalfDist> get_imag() { return dfdata_[1]; }
    std::pair<const int, const int> coord() { return coord_; }
    std::pair<const int, const int> basis() { return basis_; }

    // form_2index_small can be used with small small quantities and small large quantities
    // it must take a smal dfdist as its argument
    std::array<std::shared_ptr<Matrix>, 2> form_2index(std::shared_ptr<DFHalfComplex>);
    //std::array<std::shared_ptr<Matrix>, 2> form_2index_small(std::shared_ptr<DFHalfComplex>);
    //std::array<std::shared_ptr<Matrix>, 2> form_2index_large_large(std::shared_ptr<DFHalfComplex>);
    //std::array<std::shared_ptr<Matrix>, 2> compute_large_Jop(std::shared_ptr<const DFDist>, std::shared_ptr<const Matrix>, std::shared_ptr<const Matrix>);
    std::array<std::shared_ptr<Matrix>, 2> compute_Jop(std::shared_ptr<const DFDist>, std::shared_ptr<const Matrix>,
               std::shared_ptr<const Matrix>, std::pair<const int, const int>, std::pair<const int, const int>);
    //std::array<std::shared_ptr<Matrix>, 2> compute_small_Jop(std::shared_ptr<const DFDist>, std::array<std::shared_ptr<const Matrix>, 4>,
    //           std::array<std::shared_ptr<const Matrix>, 4>, std::pair<const int, const int>, std::pair<const int, const int>);
    std::pair<const int, const double> compute_coeff(std::pair<const int, const int>, std::pair<const int, const int>);

//    std::shared_ptr<Reference> conv_to_ref() const override;

};

}

#endif
