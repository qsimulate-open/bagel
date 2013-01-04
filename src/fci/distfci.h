//
// BAGEL - Parallel electron correlation program.
// Filename: distfci.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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

#ifndef __NEWINT_FCI_DISTFCI_H
#define __NEWINT_FCI_DISTFCI_H

#include <src/fci/harrison.h>

namespace bagel {

// Parallel FCI based on Harrison-Zarrabian algorithm.
// The algorithm is basically the same as written in    
// Z. Gan and R. J. Harrison, SC '05: Proceedings of the 2005 ACM/IEEE conference on Supercomputing
//
// The implementation is based on the HarrisonZarrabian class written by Shane Parker.

class DistFCI : public HarrisonZarrabian {

  protected:
    // const denom (so far implemented in HarrisonZarrabian)
#if 0
    void const_denom() override;
#endif

    std::shared_ptr<Dvec> form_sigma(std::shared_ptr<const Dvec> c, std::shared_ptr<const MOFile> jop, const std::vector<int>& conv) const override;

    // run-time functions implemented in HarrisonZarrabian
    void sigma_1    (std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, std::shared_ptr<const MOFile> jop) const override;
    void sigma_2bb  (std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, std::shared_ptr<const MOFile> jop) const override;
    void sigma_3    (std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, std::shared_ptr<const MOFile> jop) const override;
#if 0
    void sigma_2aa  (std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, std::shared_ptr<const MOFile> jop) const;
#endif
    // this function substitute sigma_2ab_1/2/3 in parallel execution. Dvec cannot be used here.
    void sigma_2ab(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, std::shared_ptr<const MOFile> jop) const;

  public:
    // this constructor is ugly... to be fixed some day...
    DistFCI(const std::multimap<std::string, std::string> a, std::shared_ptr<const Reference> b,
            const int ncore = -1, const int nocc = -1, const int nstate = -1);
    
#if 0
    void update(std::shared_ptr<const Coeff>) override;
#endif
};

}

#endif

