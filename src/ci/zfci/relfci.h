//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relfci.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_CI_ZFCI_RELFCI_H
#define __SRC_CI_ZFCI_RELFCI_H

#include <src/ci/zfci/zharrison.h>

namespace bagel {

class RelFCI : public ZHarrison {
  protected:
    // breit and gaunt
    bool gaunt_;
    bool breit_;

    void dump_integrals_and_exit() const override;
    std::shared_ptr<const ZCoeff_Block> init_coeff() override;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<ZHarrison>(*this) & gaunt_ & breit_;
    }

  public:
    RelFCI(std::shared_ptr<const PTree> a, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> b,
           const int ncore = -1, const int nocc = -1, std::shared_ptr<const ZCoeff_Block> coeff_zcas = nullptr, const bool store_c = false, const bool store_g = false);
    void update(std::shared_ptr<const ZCoeff_Block> coeff) override;
};

}

#endif
