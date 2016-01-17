//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/dist_form_sigma.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __BAGEL_RAS_DIST_FORM_SIGMA_RAS_H
#define __BAGEL_RAS_DIST_FORM_SIGMA_RAS_H

#include <src/ci/ras/civector.h>
#include <src/ci/fci/mofile.h>

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

