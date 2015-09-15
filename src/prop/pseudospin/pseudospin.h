//
// BAGEL - Parallel electron correlation program.
// Filename: pseudospin.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynolds2018@u.northwestern.edu>
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


#ifndef __SRC_PROP_PSEUDOSPIN_PSEUDOSPIN_H
#define __SRC_PROP_PSEUDOSPIN_PSEUDOSPIN_H

#include <src/util/math/zmatrix.h>

namespace bagel {

// Some operator that acts on a basis of pseudospin eigenstates
class Spin_Operator {
  protected:
    int nspin_;
    std::shared_ptr<ZMatrix> matrix_;
    std::string coeff_name_;

  public:
    Spin_Operator(const std::shared_ptr<ZMatrix> mat, const std::string name)
     : matrix_(mat), coeff_name_(name) {
      assert(matrix_->ndim() == matrix_->mdim());
      nspin_ = matrix_->ndim();
    }

    int nspin() const { return nspin_; }
    std::shared_ptr<ZMatrix> matrix() const { return matrix_; }
    std::string coeff_name() const { return coeff_name_; }
};

class Pseudospin {
  protected:

  public:
    void compute_extended_stevens_operators() const;

};

}

#endif
