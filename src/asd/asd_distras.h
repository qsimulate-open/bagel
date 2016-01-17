//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_distras.h
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __ASD_ASD_DISTRAS_H
#define __ASD_ASD_DISTRAS_H

#include <src/asd/asd.h>

namespace bagel {

class ASD_DistRAS : public ASD<DistRASDvec> {
  public:
    ASD_DistRAS(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DimerDistRAS> cispace)
      : ASD<DistRASDvec>(input, dimer, cispace) {}

  private:
    std::shared_ptr<DistRASDvec> form_sigma(std::shared_ptr<const DistRASDvec> ccvec, std::shared_ptr<const MOFile> jop) const override;
    std::shared_ptr<DistRASDvec> form_sigma_1e(std::shared_ptr<const DistRASDvec> ccvec, const double* modata) const override;

    void sigma_aa(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASCivec> sigma, const double* h1, const double* h2) const;
    void sigma_bb(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASCivec> sigma, const double* h1, const double* h2) const;
    void sigma_ab(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASCivec> sigma, const double* h1, const double* h2) const;

    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_rdm12_monomer(std::shared_ptr<const DistRASDvec> civec, const int i) const override;

    std::shared_ptr<DistRASDvec> contract_I(std::shared_ptr<const DistRASDvec> A, std::shared_ptr<Matrix> coef, int offset, int nstA, int nstB, int nstates) const override;
    std::shared_ptr<DistRASDvec> contract_J(std::shared_ptr<const DistRASDvec> B, std::shared_ptr<Matrix> coef, int offset, int nstA, int nstB, int nstates) const override;
};

}

#endif
