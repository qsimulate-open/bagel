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

#ifndef __BAGEL_FCI_DISTFCI_H
#define __BAGEL_FCI_DISTFCI_H

#include <memory>
#include <vector>
#include <src/fci/fci.h>
#include <src/fci/space.h>

namespace bagel {

// Parallel FCI based on Harrison-Zarrabian algorithm.
// The algorithm is basically the same as written in
// Z. Gan and R. J. Harrison, SC '05: Proceedings of the 2005 ACM/IEEE conference on Supercomputing
//
// The implementation is based on the HarrisonZarrabian class written by Shane Parker.

class DistFCI : public FCI {

  protected:
    std::shared_ptr<Space> space_;
    std::shared_ptr<DistCivec> denom_;

    // const_denom function here only makes a denom for local data of DistCivec. This is called from FCI
    void const_denom() override;

    // denominator
    void generate_guess(const int nspin, const int nstate, std::vector<std::shared_ptr<DistCivec>> out);

    // Determinant seeds in parallel
    std::vector<std::pair<std::bitset<nbit__>, std::bitset<nbit__>>> detseeds(const int ndet) override;

    // just to confort FCI
    std::shared_ptr<Dvec> form_sigma(std::shared_ptr<const Dvec> c, std::shared_ptr<const MOFile> jop, const std::vector<int>& conv) const override
      { throw std::logic_error("DistFCI::form_sigma disabled"); }

    std::vector<std::shared_ptr<DistCivec>> form_sigma(std::vector<std::shared_ptr<DistCivec>>&, std::shared_ptr<const MOFile> jop, const std::vector<int>& conv) const;

    void sigma_bb(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop,
                  const std::shared_ptr<const Determinants>, const std::shared_ptr<const Determinants>) const;

    void sigma_bb(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop) const;
    void sigma_aa(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop) const;
    void sigma_ab(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop) const;

  public:
    // this constructor is ugly... to be fixed some day...
    DistFCI(const std::shared_ptr<const PTree> a, std::shared_ptr<const Reference> b,
            const int ncore = -1, const int nocc = -1, const int nstate = -1);

    void update(std::shared_ptr<const Coeff>) override;

    // FCI compute function using DistCivec
    void compute() override;

};

}

#endif

