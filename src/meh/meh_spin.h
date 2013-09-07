//
// BAGEL - Parallel electron correlation program.
// Filename: meh_spin.h
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifndef __meh_meh_spin_h
#define __meh_meh_spin_h

#include <src/math/sparsematrix.h>

namespace bagel {

class MEHSpin : public SparseMatrix {
  protected:
    const int max_spin_;

  public:
    MEHSpin(const int dimension, std::map<std::pair<int,int>, double>& coords, const int max) : SparseMatrix(dimension, dimension, coords), max_spin_(max) {}

    const int max() const { return max_spin_; }

    std::shared_ptr<Matrix> apply(const Matrix& o) const;
    void filter(Matrix& o, const int desired_spin) const;
};

}

#endif
