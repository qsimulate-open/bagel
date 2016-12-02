//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhcore.h
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Raymond Wang <yiqunwang2021@u.northwestern.edu> 
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


#ifndef __SRC_DKH_DKHCORE_H
#define __SRC_DKH_DKHCORE_H

#include <src/util/constants.h>
#include <src/util/math/matrix.h>
#include <src/util/math/zmatrix.h>
#include <src/molecule/molecule.h>

namespace bagel {

template<typename DataType>
class DKHcore_ : public std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type {
  protected:
    using MatType = typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type;
    std::shared_ptr<const Molecule> mol_;

    void init(std::shared_ptr<const Molecule>) { assert(false); }
  public:
    DKHcore_() { }
    DKHcore_(std::shared_ptr<const Molecule>);
};

template<>
void DKHcore_<double>::init(std::shared_ptr<const Molecule>);

using DKHcore = DKHcore_<double>;
using ZDKHcore = DKHcore_<std::complex<double>>;

extern template class DKHcore_<double>;
extern template class DKHcore_<std::complex<double>>;

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::DKHcore_<double>)
BOOST_CLASS_EXPORT_KEY(bagel::DKHcore_<std::complex<double>>)

#endif

