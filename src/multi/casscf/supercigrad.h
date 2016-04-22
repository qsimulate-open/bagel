//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: supercigrad.h
// Copyright (C) 2011 Toru Shiozaki
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


#ifndef __BAGEL_CASSCF_SUPERCIGRAD_H
#define __BAGEL_CASSCF_SUPERCIGRAD_H

#include <src/multi/casscf/superci.h>

namespace bagel {

class SuperCIGrad : public SuperCI {
  public:
    SuperCIGrad(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref)
      : SuperCI(idat, geom, ref) { }

    void compute() {
      // compute CASSCF fist
      SuperCI::compute();

      // then make sure that the orbitals are natural orbitals
      // in the worst case, coeff and RDMs are not consistent and coeff is not natural orbital...
      fci_->compute();
      fci_->compute_rdm12();

      // form natural orbitals and update coefficients
      const std::pair<std::shared_ptr<Matrix>, VectorB> natorb = fci_->natorb_convert();
      coeff_ = update_coeff(coeff_, natorb.first);
      occup_ = natorb.second;
      if (natocc_) print_natocc();
    }

};

}

#endif
