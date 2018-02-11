//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hcoreinfo.h
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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


#ifndef __SRC_WFN_HCOREINFO_H
#define __SRC_WFN_HCOREINFO_H

#include <src/molecule/molecule.h>

// Contains info on Hcore.
//
//   TODO
//      o derived class for each method (DKH, ECP, ...)
//      o for the time being, ECP parameters are saved in Atom object.
//        Ideally, they should be placed here.

namespace bagel {

enum HcoreType { standard, dkh, ecp };

class HcoreInfo {
  protected:
    HcoreType type_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & type_;
    }

  public:
    HcoreInfo() : type_(HcoreType::standard) { }
    HcoreInfo(std::shared_ptr<const PTree> idata);

    bool dkh() const { return type_ == HcoreType::dkh; }
    bool ecp() const { return type_ == HcoreType::ecp; }
    bool standard() const { return type_ == HcoreType::standard; }
    void print() const;

    // DKH specific
    std::shared_ptr<Matrix> compute_dkh(std::shared_ptr<const Molecule> current) const;

    // general function
    std::shared_ptr<Matrix> compute(std::shared_ptr<const Molecule> current) const;
};

}

#endif
