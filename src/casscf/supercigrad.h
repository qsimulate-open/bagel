//
// Newint - Parallel electron correlation program.
// Filename: supercigrad.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __NEWINT_CASSCF_SUPERCIGRAD_H
#define __NEWINT_CASSCF_SUPERCIGRAD_H

#include <memory>
#include <string>
#include <map>
#include <src/scf/scf.h>
#include <src/casscf/superci.h>
#include <src/wfn/rdm.h>

class SuperCIGrad : public SuperCI {

  protected:

  public:
    SuperCIGrad(const std::multimap<std::string, std::string> idat, const std::shared_ptr<const Geometry> geom)
      : SuperCI(idat, geom) { };
    ~SuperCIGrad() {};

    void compute() {
      SuperCI::compute();
      fci_->compute_rdm12();
      form_natural_orbs();
      fci_->update(coeff_);
      fci_->compute();
      fci_->compute_rdm12();
    };

};

#endif
