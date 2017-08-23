//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fock_base.h
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


#ifndef __SRC_SCF_FOCK_BASE_H
#define __SRC_SCF_FOCK_BASE_H

#include <src/wfn/geometry.h>
#include <src/mat1e/matrix1e.h>

namespace bagel {

template <typename MatType = Matrix, class Enable = typename std::enable_if<(std::is_same<MatType, Matrix>::value || std::is_same<MatType, ZMatrix>::value)>::type>
class Fock_base_ : public Matrix1e_<MatType, Enable> {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const MatType> previous_;
    std::shared_ptr<const MatType> density_;
    void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int, std::shared_ptr<const Molecule>) override;

    // virtual function that is to be defined in the derived class
    virtual void fock_two_electron_part(std::shared_ptr<const MatType>) = 0;
    void fock_one_electron_part();

    // for non-DF Fock builds
    std::vector<double> schwarz_;
    double schwarz_thresh_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Matrix1e_<MatType, Enable>>(*this) & geom_ & previous_ & density_ & schwarz_ & schwarz_thresh_;
    }

  public:
    Fock_base_() { }
    Fock_base_(std::shared_ptr<const Geometry>, std::shared_ptr<const MatType> prev, std::shared_ptr<const MatType> den,
               const std::vector<double>& = std::vector<double>());
    virtual ~Fock_base_() { }

    using MatType::mdim;
    using MatType::ndim;
    using MatType::element;
    using MatType::fill_upper_conjg;
    using Matrix1e_<MatType, Enable>::init;

};

using Fock_base = Fock_base_<Matrix>;
using Fock_base_London = Fock_base_<ZMatrix>;

}

extern template class bagel::Fock_base_<bagel::Matrix>;
extern template class bagel::Fock_base_<bagel::ZMatrix>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Fock_base)
BOOST_CLASS_EXPORT_KEY(bagel::Fock_base_London)

#endif
