//
// BAGEL - Parallel electron correlation program.
// Filename: ras/dist_form_sigma.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __BAGEL_RAS_DIST_FORM_SIGMA_RAS_H
#define __BAGEL_RAS_DIST_FORM_SIGMA_RAS_H

#include <src/ras/civector.h>
#include <src/fci/mofile.h>

namespace bagel {

class DistFormSigmaRAS {
  protected:
    const bool sparse_;

  public:
    DistFormSigmaRAS(const bool sparse) : sparse_(sparse) {}

    // This is really all this class is
    std::shared_ptr<DistRASDvec> operator()(std::shared_ptr<const DistRASDvec> cc, std::shared_ptr<const MOFile> jop, const std::vector<int>& conv) const;
    std::shared_ptr<DistRASDvec> operator()(std::shared_ptr<const DistRASDvec> cc, const double* mo1e) const;

  private:
    // Helper functions for sigma formation
    void sigma_bb(std::shared_ptr<const DistRASCivec> cc, std::shared_ptr<DistRASCivec> sigma, const double* g, const double* mo2e) const;
    void sigma_ab(std::shared_ptr<const DistRASCivec> cc, std::shared_ptr<DistRASCivec> sigma, const double* mo2e) const;
};

}

#endif

