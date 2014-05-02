//
// BAGEL - Parallel electron correlation program.
// Filename: zcoeff.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __src_london_zcoeff_h
#define __src_london_zcoeff_h

#include <src/wfn/geometry_london.h>
#include <src/math/zmatrix.h>

namespace bagel {

class ZCoeff : public ZMatrix {
  protected:
    std::shared_ptr<const Geometry_London> cgeom_;

  private:
    int num_basis(std::vector<std::shared_ptr<const ZCoeff>> coeff_vec) const;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<ZMatrix>(*this) & cgeom_;
    }

  public:
    ZCoeff() { }
    ZCoeff(const ZMatrix&);
    ZCoeff(std::vector<std::shared_ptr<const ZCoeff>> coeff_vec);
    ZCoeff(std::shared_ptr<const Geometry_London> g) : ZMatrix(g->nbasis(), g->nbasis()), cgeom_(g) {}

    std::shared_ptr<const Geometry_London> geom() const { assert(cgeom_); return cgeom_; }

    std::shared_ptr<ZMatrix> form_weighted_density_rhf(const int n, const std::vector<std::complex<double>>& e, const int offset = 0) const;
    std::pair<std::shared_ptr<ZMatrix>, std::shared_ptr<ZMatrix>> split(const int, const int) const;
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::ZCoeff)

#endif
