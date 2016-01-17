//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: soscf.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __BAGEL_SRC_SCF_SOSCF_H
#define __BAGEL_SRC_SCF_SOSCF_H

#include <src/scf/scf_base.h>
#include <src/mat1e/sohcore.h>

namespace bagel {

class SOSCF : public SCF_base {
  protected:

    bool dodf_;

    std::shared_ptr<const SOHcore> sohcore_;
    std::shared_ptr<const ZMatrix> socoeff_;
    std::shared_ptr<const ZMatrix> sooverlap_;
    std::shared_ptr<const ZMatrix> sotildex_;
    VectorB soeig_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<SCF_base>(*this);
      ar & dodf_ & sohcore_ & socoeff_ & sooverlap_ & socoeff_ & sotildex_;
    }

  public:
    SOSCF() { }
    SOSCF(const std::shared_ptr<const PTree> idata_, const std::shared_ptr<const Geometry> geom,
          const std::shared_ptr<const Reference> re = nullptr);

    void compute() override;

    std::shared_ptr<const ZMatrix> sooverlap();
    std::shared_ptr<const ZMatrix> sotildex();
    VectorB& soeig() { return soeig_; }

    bool dodf() const { return dodf_; }

    std::shared_ptr<const Reference> conv_to_ref() const override;
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::SOSCF)

#endif
