//
// BAGEL - Parallel electron correlation program.
// Filename: dmp2.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __SRC_REL_RELMP2_H
#define __SRC_REL_RELMP2_H

#include <src/rel/dirac.h>
#include <src/wfn/method.h>

namespace bagel {

class DMP2 : public Method {
  protected:
    std::shared_ptr<Dirac> scf_;
    int ncore_;

    std::string abasis_;

    double energy_;

  public:
    DMP2(const std::shared_ptr<const PTree>, const std::shared_ptr<const Geometry>,
         const std::shared_ptr<const Reference> = std::shared_ptr<const Reference>());

    virtual void compute() override;
    virtual std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; } 

    double energy() const { return energy_; }
    int ncore() const { return ncore_; }
    std::string abasis() const { return abasis_; }
};

}

#endif
