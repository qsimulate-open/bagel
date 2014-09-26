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

template <typename T, class U> class Matrix1eTask_;

// specialized matrix for 1e integrals
template <typename MatType = Matrix, class Enable = typename std::enable_if<(std::is_same<MatType, Matrix>::value || std::is_same<MatType, ZMatrix>::value)>::type>
class Matrix1e_ : public MatType{

  friend class Matrix1eTask_<MatType, Enable>;
  protected:
    virtual void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int, std::shared_ptr<const Molecule>) = 0;
    virtual void init(std::shared_ptr<const Molecule>);

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<MatType>(*this);
    }

  public:
    Matrix1e_() { }
    Matrix1e_(const std::shared_ptr<const Molecule>);
    Matrix1e_(const Matrix1e_&);
    virtual ~Matrix1e_() { }

    using MatType::zero;
    using MatType::size;
    using MatType::data;
    using MatType::allreduce;

};

using Matrix1e = Matrix1e_<Matrix>;
using ZMatrix1e = Matrix1e_<ZMatrix>;

}

extern template class bagel::Matrix1e_<bagel::Matrix>;
extern template class bagel::Matrix1e_<bagel::ZMatrix>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Matrix1e)
BOOST_CLASS_EXPORT_KEY(bagel::ZMatrix1e)

#endif
