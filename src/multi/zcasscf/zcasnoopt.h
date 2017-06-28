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
#include <src/prop/pseudospin/pseudospin.h>

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
      if (external_rdm_.empty()) {
        std::cout << " Computing RDMs from FCI calculation " << std::endl;
        fci_->compute();
        fci_->compute_rdm12();
      } else {
        fci_->read_external_rdm12_av(external_rdm_);
      }
      energy_ = fci_->energy();

      // TODO When the Property class is implemented, this should be one
      std::shared_ptr<const PTree> aniso_data = idata_->get_child_optional("aniso");
      if (aniso_data) {
        if (geom_->magnetism()) {
          std::cout << "  ** Magnetic anisotropy analysis is currently only available for zero-field calculations; sorry." << std::endl;
        } else {
          const int nspin = aniso_data->get<int>("nspin", (idata_->get_vector<int>("state", 0)).size()-1);
          Pseudospin ps(nspin, geom_, fci_->conv_to_ciwfn(), aniso_data);
          ps.compute(energy_, coeff_->active_part());
        }
      }
    }

};

}

#endif
