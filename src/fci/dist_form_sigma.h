//
// BAGEL - Parallel electron correlation program.
// Filename: fci/dist_form_sigma.h
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

#ifndef __BAGEL_FCI_DIST_FORM_SIGMA_H
#define __BAGEL_FCI_DIST_FORM_SIGMA_H

#include <src/fci/civec.h>
#include <src/fci/mofile.h>
#include <src/fci/space.h>

namespace bagel {

class FormSigmaDistFCI {
  protected:
    std::shared_ptr<const Space_base> space_;

  public:
    FormSigmaDistFCI(std::shared_ptr<const Space_base> sp = std::shared_ptr<const Space_base>()) : space_(sp) {}

    std::vector<std::shared_ptr<DistCivec>> operator()(const std::vector<std::shared_ptr<DistCivec>>& cc, std::shared_ptr<const MOFile> jop, const std::vector<int>& conv) const;
    std::shared_ptr<DistDvec> operator()(std::shared_ptr<const DistDvec> cc, std::shared_ptr<const MOFile> jop) const;
    //std::shared_ptr<DistDvec> operator()(std::shared_ptr<const DistDvec> cc, const double* mo1e) const;

  private:
    // Helper functions for sigma formation
    void sigma_bb(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop,
                  const std::shared_ptr<const Determinants>, const std::shared_ptr<const Determinants>) const;

    void sigma_bb(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop) const;
    void sigma_aa(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop, std::shared_ptr<const Determinants> int_det) const;
    void sigma_ab(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop) const;
};

}

#endif

