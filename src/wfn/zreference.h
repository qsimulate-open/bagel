//
// BAGEL - Parallel electron correlation program.
// Filename: zreference.h
// Copyright (C) 2015 Toru Shiozaki
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
    ZReference(std::shared_ptr<const Geometry> g, std::shared_ptr<const ZCoeff> c, const double en, const int nocc, const int nact, const int nvirt)
     : Reference(g, nullptr, nocc, nact, nvirt, en), zcoeff_(c) {
      assert(geom_->magnetism());  // Currently only used for GIAO wavefunctions
    }

    std::shared_ptr<const Coeff> coeff() const override { throw std::logic_error("ZReference::coeff() should not be called"); }
    std::shared_ptr<const ZCoeff> zcoeff() const { return zcoeff_; }

    std::shared_ptr<Reference> project_coeff(std::shared_ptr<const Geometry> geomin, const bool check_geom_change = true) const override;

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::ZReference)

#endif
