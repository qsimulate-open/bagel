//
// BAGEL - Parallel electron correlation program.
// Filename: supercigrad.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __BAGEL_CASSCF_SUPERCIGRAD_H
#define __BAGEL_CASSCF_SUPERCIGRAD_H

#include <src/casscf/superci.h>

namespace bagel {

class SuperCIGrad : public SuperCI {

  protected:

  public:
    SuperCIGrad(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref)
      : SuperCI(idat, geom, ref) { };

    void compute() {
      // compute CASSCF fist
      SuperCI::compute();

      // then make sure that the orbitals are natural orbitals
      // in the worst case, coeff and RDMs are not consistent and coeff is not natural orbital...
      fci_->update(coeff_);
      fci_->compute();
      fci_->compute_rdm12();
// TODO form_naturla_orbs might have a bug...
      form_natural_orbs();
      fci_->update(coeff_);
      fci_->compute();
      fci_->compute_rdm12();
      form_natural_orbs();
      fci_->update(coeff_);
      fci_->compute();
      fci_->compute_rdm12();
    };

};

}

#endif
