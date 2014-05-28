//
// BAGEL - Parallel electron correlation program.
// Filename: meh_distcas.h
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __MEH_MEH_DIST_CAS_H
#define __MEH_MEH_DIST_CAS_H

#include <src/meh/meh.h>
#include <src/fci/civec.h>

namespace bagel {

class MEH_DistCAS : public MultiExcitonHamiltonian<DistDvec> {
   public:
      MEH_DistCAS(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DimerDistCAS> cispace) :
        MultiExcitonHamiltonian<DistDvec>(input, dimer, cispace)
      {
        cispace_->intermediates();
      }

   private:
      std::shared_ptr<DistDvec> form_sigma(std::shared_ptr<const DistDvec> ccvec, std::shared_ptr<const MOFile> jop) const override;
      std::shared_ptr<DistDvec> form_sigma_1e(std::shared_ptr<const DistDvec> ccvec, const double* modata) const override;
};

}

#endif
