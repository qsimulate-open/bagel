//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/form_sigma.h
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

#ifndef __BAGEL_RAS_FORM_SIGMA_RAS_H
#define __BAGEL_RAS_FORM_SIGMA_RAS_H

#include <src/ci/ras/civector.h>
#include <src/ci/fci/mofile.h>

namespace bagel {

class FormSigmaRAS {
  protected:
    int batchsize_;

  public:
    FormSigmaRAS(const int b = 512) : batchsize_(b) {}

    /// Applies Hamiltonian to cc using the provided MOFile, skipping the vectors marked as converged
    std::shared_ptr<RASDvec> operator()(std::shared_ptr<const RASDvec> ccvec, std::shared_ptr<const MOFile> jop, const std::vector<int>& conv) const;
    /// Applies Hamiltonian to cc using the provided 1e and 2e integrals, skipping the vectors marked as converged
    std::shared_ptr<RASDvec> operator()(std::shared_ptr<const RASDvec> ccvec, std::shared_ptr<const Matrix> mo1e, std::shared_ptr<const Matrix> mo2e, const std::vector<int>& conv) const;
    void operator()(const RASCivecView cc, RASCivecView sigma, std::shared_ptr<const MOFile> jop) const;

  private:
    // Helper functions for sigma formation
    void sigma_aa(const RASCivecView cc, RASCivecView sigma, const double* g, const double* mo2e) const;
    void sigma_bb(const RASCivecView cc, RASCivecView sigma, const double* g, const double* mo2e) const;
    void sigma_ab(const RASCivecView cc, RASCivecView sigma, const double* mo2e) const;
};

}

#endif

