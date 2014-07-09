//
// BAGEL - Parallel electron correlation program.
// Filename: relreference_london.h
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

#ifndef __SRC_LONDON_RELREFERENCE_LONDON_H
#define __SRC_LONDON_RELREFERENCE_LONDON_H

#include <src/london/reference_london.h>

namespace bagel {

class RelReference_London : public Reference_London {
  protected:
    bool gaunt_;
    bool breit_;
    int nneg_;
    std::shared_ptr<const ZMatrix> relcoeff_;
    std::shared_ptr<const ZMatrix> relcoeff_full_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Reference_London>(*this) & gaunt_ & breit_ & nneg_ & relcoeff_ & relcoeff_full_;
    }

  public:
    RelReference_London() { }
    RelReference_London(std::shared_ptr<const Geometry_London> g, std::shared_ptr<const ZMatrix> c, const double en, const int nneg, const int nocc,
                        const int nvirt, const bool ga, const bool br)
     : Reference_London(g, nullptr, nocc, 0, nvirt, en), gaunt_(ga), breit_(br), nneg_(nneg), relcoeff_(c->slice_copy(nneg_, c->mdim())), relcoeff_full_(c) {
    }

    const std::shared_ptr<const ZCoeff> zcoeff() const override { throw std::logic_error("RelReference_London::zcoeff() should not be called"); }
    const std::shared_ptr<const ZMatrix> relcoeff() const { return relcoeff_; }
    const std::shared_ptr<const ZMatrix> relcoeff_full() const { return relcoeff_full_; }

    bool gaunt() const { return gaunt_; }
    bool breit() const { return breit_; }
    int nneg() const { return nneg_; }

    std::shared_ptr<Reference> project_coeff(std::shared_ptr<const Geometry_London> geomin) const override;

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RelReference_London)

#endif
