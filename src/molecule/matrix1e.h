//
// BAGEL - Parallel electron correlation program.
// Filename: matrix1e.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#ifndef __SRC_MOLECULE_MATRIX1E_H
#define __SRC_MOLECULE_MATRIX1E_H

#include <src/molecule/molecule.h>
#include <src/math/matrix.h>

namespace bagel {

// specialized matrix for 1e integrals
class Matrix1e : public Matrix {
  protected:
    std::shared_ptr<const Molecule> mol_;

    virtual void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int) = 0;
    virtual void init();

  public:
    Matrix1e(const std::shared_ptr<const Molecule>);
    Matrix1e(const std::shared_ptr<const Molecule>, const int n, const int m);
    Matrix1e(const Matrix1e&);

    const std::shared_ptr<const Molecule> mol() const { return mol_; };

};

}

#endif
