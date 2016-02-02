//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci/dist_form_sigma.h
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

#ifndef __BAGEL_FCI_DIST_FORM_SIGMA_H
#define __BAGEL_FCI_DIST_FORM_SIGMA_H

#include <src/ci/fci/distcivec.h>
#include <src/ci/fci/mofile.h>
#include <src/ci/fci/space.h>

namespace bagel {

class FormSigmaDistFCI {
  protected:
    std::shared_ptr<const Space_base> space_;

  public:
    FormSigmaDistFCI(std::shared_ptr<const Space_base> sp = nullptr) : space_(sp) {}

    std::vector<std::shared_ptr<DistCivec>> operator()(const std::vector<std::shared_ptr<DistCivec>>& cc, std::shared_ptr<const MOFile> jop, const std::vector<int>& conv) const;
    std::shared_ptr<DistDvec> operator()(std::shared_ptr<const DistDvec> cc, std::shared_ptr<const MOFile> jop) const;
    //std::shared_ptr<DistDvec> operator()(std::shared_ptr<const DistDvec> cc, const double* mo1e) const;

  private:
    // Helper functions for sigma formation
    void sigma_bb(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop,
                  std::shared_ptr<const Determinants>, std::shared_ptr<const Determinants>) const;

    void sigma_bb(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop) const;
    void sigma_aa(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop, std::shared_ptr<const Determinants> int_det) const;
    void sigma_ab(std::shared_ptr<const DistCivec> cc, std::shared_ptr<DistCivec> sigma, std::shared_ptr<const MOFile> jop) const;
};

}

#endif

