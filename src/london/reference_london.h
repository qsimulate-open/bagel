//
// BAGEL - Parallel electron correlation program.
// Filename: reference_london.h
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

#ifndef __SRC_LONDON_REFERENCE_LONDON_H
#define __SRC_LONDON_REFERENCE_LONDON_H

#include <src/wfn/reference.h>
#include <src/molecule/zhcore.h>

namespace bagel {

class Reference_London : public Reference {
  protected:
    std::shared_ptr<const ZCoeff> zcoeff_;
    std::shared_ptr<const ZHcore> zhcore_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Reference>(*this) & zcoeff_ & zhcore_;
    }

  public:
    Reference_London() { }
    Reference_London(std::shared_ptr<const Geometry> g, std::shared_ptr<const ZCoeff> c,
                     const int nclo, const int nact, const int nvirt, const double en = 0.0);

    std::shared_ptr<const Coeff> coeff() const override { throw std::logic_error("Reference_London::coeff() should not be called"); }
    std::shared_ptr<const Hcore> hcore() const override { throw std::logic_error("Reference_London::hcore() should not be called"); }
    std::shared_ptr<Reference> project_coeff(std::shared_ptr<const Geometry> geomin) const override;

    virtual std::shared_ptr<const ZCoeff> zcoeff() const { return zcoeff_; }

    std::shared_ptr<const ZHcore> zhcore() const { return zhcore_; }

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Reference_London)

#endif
