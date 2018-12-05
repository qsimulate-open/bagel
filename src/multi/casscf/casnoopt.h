//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: casnoopt.h
// Copyright (C) 2016 Toru Shiozaki
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

#ifndef __SRC_MULTI_CASSCF_CASNOOPT_H
#define __SRC_MULTI_CASSCF_CASNOOPT_H

#include <src/multi/casscf/casscf.h>
#include <src/prop/hyperfine.h>

namespace bagel {

class CASNoopt : public CASSCF {
  protected:
    void common_init() {
      std::cout << "    * No orbital optimization will be performed!" << std::endl << std::endl;
    }

  public:
    CASNoopt(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref)
     : CASSCF(idat, geom, ref) { common_init(); }

    void compute() override {
      if (nact_) {
        if (external_rdm_.empty()) {
          fci_->compute();
          fci_->compute_rdm12();
        } else {
          fci_->read_external_rdm12_av(external_rdm_);
        }
        energy_ = fci_->energy();
      } else {
        energy_ = ref_->energy();
      }

      if (canonical_) {
        auto tmp = semi_canonical_orb();
        coeff_ = std::get<0>(tmp);
        eig_   = std::get<1>(tmp);
        occup_ = std::get<2>(tmp);
      }

      if (do_hyperfine_ && !geom_->external() && nstate_ == 1 && external_rdm_.empty()) {
        HyperFine hfcc(geom_, spin_density(), fci_->det()->nspin(), "CASSCF");
        hfcc.compute();
      }
    }

    std::shared_ptr<const Reference> conv_to_ref() const override {
      std::shared_ptr<Reference> out;
      if (nact_) {
        out = std::make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, energy_,
                                          fci_->rdm1(), fci_->rdm2(), fci_->rdm1_av(), fci_->rdm2_av(), fci_->conv_to_ciwfn());
      } else {
        out = std::make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, energy_);
      }
      out->set_eig(eig_);
      out->set_occup(occup_);
      return out;
    }
};

}

#endif
