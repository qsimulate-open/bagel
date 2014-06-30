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


#ifndef __SRC_MOLECULE_MATRIX1E_H
#define __SRC_MOLECULE_MATRIX1E_H

#include <src/molecule/molecule.h>
#include <src/math/matrix.h>

namespace bagel {

// specialized matrix for 1e integrals
class Matrix1e : public Matrix {
  friend class Matrix1eTask;
  protected:
    std::shared_ptr<Matrix> soab_, soba_, soaa_;
    virtual void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int, std::shared_ptr<const Molecule>) = 0;
    virtual void init(std::shared_ptr<const Molecule>);

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Matrix>(*this);
    }

  public:
    Matrix1e() { }
    Matrix1e(const std::shared_ptr<const Molecule>);
    Matrix1e(const Matrix1e&);
    virtual ~Matrix1e() { }

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Matrix1e)

#endif
