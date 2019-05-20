//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: matrix1e.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_MAT1E_MATRIX1E_H
#define __SRC_MAT1E_MATRIX1E_H

#include <src/molecule/molecule.h>
#include <src/util/math/matrix.h>

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
    Matrix1e_(std::shared_ptr<const Molecule>);
    Matrix1e_(const Matrix1e_&);
    virtual ~Matrix1e_() { }

    Matrix1e_& operator=(const Matrix1e_&) = default;
    Matrix1e_& operator=(Matrix1e_&&) = default;

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
