//
// BAGEL - Parallel electron correlation program.
// Filename: dimer_scf.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Modified by: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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


#ifndef __BAGEL_DIMER_DIMER_SCF_H
#define __BAGEL_DIMER_DIMER_SCF_H

#include <src/scf/scf_base.h>
#include <src/dimer/dimer_levelshift.h>
#include <src/dimer/dimer.h>

namespace bagel {

class DimerSCF : public SCF_base {
  protected:
    std::shared_ptr<const Dimer> dimer_;
    std::shared_ptr<ShiftDimer> levelshift_;

  public:
    DimerSCF(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Dimer> dimer);

    void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override;
};

}

#endif
