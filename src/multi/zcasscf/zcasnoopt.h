//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcasnoopt.h
// Copyright (C) 2016 Toru Shiozaki
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

#ifndef __SRC_MULTI_ZCASSCF_ZCASNOOPT_H
#define __SRC_MULTI_ZCASSCF_ZCASNOOPT_H

#include <src/multi/zcasscf/zcasscf.h>

namespace bagel {

class ZCASNoopt : public ZCASSCF {
  protected:
    void common_init() {
      std::cout << "    * No orbital optimization will be performed!" << std::endl << std::endl;
    }

  public:
    ZCASNoopt(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref)
     : ZCASSCF(idat, geom, ref) { common_init(); }

    void compute() override {
      Timer fci_time(0);
      if (external_rdm_.empty()) {
        std::cout << " Computing RDMs from FCI calculation " << std::endl;
        fci_->compute();
        fci_->compute_rdm12();
        fci_time.tick_print("ZFCI");
      } else {
        fci_->read_external_rdm12_av(external_rdm_);
      }
      fci_time.tick_print("RDMs");
      energy_ = fci_->energy();
    }

};

}

#endif
