//
// BAGEL - Parallel electron correlation program.
// Filename: coeff.h
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


#ifndef __SRC_SCF_COEFF_H
#define __SRC_SCF_COEFF_H

#include <src/wfn/geometry.h>

namespace bagel {

template <typename MatType = Matrix>
class Coeff_ : public MatType {
  private:
    int num_basis(std::vector<std::shared_ptr<const Coeff_<MatType>>> coeff_vec) const;

    using DataType = typename MatType::value_type;

    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<MatType>(*this);
    }

  public:
    Coeff_() { }
    Coeff_(const MatType&);
    Coeff_(MatType&&);
    Coeff_(std::vector<std::shared_ptr<const Coeff_<MatType>>> coeff_vec);
    Coeff_(std::shared_ptr<const Geometry> g) : MatType(g->nbasis(), g->nbasis()) {}

    std::shared_ptr<MatType> form_weighted_density_rhf(const int n, const std::vector<double>& e) const;
    std::pair<std::shared_ptr<MatType>, std::shared_ptr<MatType>> split(const int, const int) const;

  public:
    using MatType::data;
    using MatType::mdim;
    using MatType::ndim;
    using MatType::size;
    using MatType::slice;
};

using Coeff = Coeff_<Matrix>;
using ZCoeff = Coeff_<ZMatrix>;

}

extern template class bagel::Coeff_<bagel::Matrix>;
extern template class bagel::Coeff_<bagel::ZMatrix>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Coeff)
BOOST_CLASS_EXPORT_KEY(bagel::ZCoeff)

#endif
