//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zreference.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

// Specialized Reference object for complex orbital coefficients
// Currently used by RHF_London and SOSCF

#ifndef __SRC_WFN_ZREFERENCE_H
#define __SRC_WFN_ZREFERENCE_H

#include <src/wfn/reference.h>

namespace bagel {

class ZReference : public Reference {
  protected:
    std::shared_ptr<const ZCoeff> zcoeff_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Reference>(*this) & zcoeff_;
    }

  public:
    ZReference() { }
    ZReference(std::shared_ptr<const Geometry> g, std::shared_ptr<const ZCoeff> c, const std::vector<double> en, const int nocc, const int nact, const int nvirt)
     : Reference(g, nullptr, nocc, nact, nvirt, en), zcoeff_(c) {
    }

    // if only given one energy
    ZReference(std::shared_ptr<const Geometry> g, std::shared_ptr<const ZCoeff> c, const double en, const int nocc, const int nact, const int nvirt) :
      ZReference(g, c, std::vector<double>(1, en), nocc, nact, nvirt) { }

    // constructing a skelton
    ZReference(std::shared_ptr<const Geometry> g) : Reference(g) { }

    std::shared_ptr<const ZCoeff> zcoeff() const { return zcoeff_; }
    void set_zcoeff(std::shared_ptr<const ZMatrix> matrix) { zcoeff_ = std::make_shared<const ZCoeff>(*matrix); }

    std::shared_ptr<Reference> project_coeff(std::shared_ptr<const Geometry> geomin, const bool check_geom_change = true) const override;

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::ZReference)

#endif
