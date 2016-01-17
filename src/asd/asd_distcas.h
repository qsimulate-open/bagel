//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_distcas.h
// Copyright (C) 2013 Shane Parker
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

#ifndef __ASD_ASD_DIST_CAS_H
#define __ASD_ASD_DIST_CAS_H

#include <src/asd/asd.h>
#include <src/ci/fci/civec.h>

namespace bagel {

class ASD_DistCAS : public ASD<DistDvec> {
  public:
    ASD_DistCAS(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DimerDistCAS> cispace)
      : ASD<DistDvec>(input, dimer, cispace) {
      cispace_->intermediates();
    }

  private:
    std::shared_ptr<DistDvec> form_sigma(std::shared_ptr<const DistDvec> ccvec, std::shared_ptr<const MOFile> jop) const override;
    std::shared_ptr<DistDvec> form_sigma_1e(std::shared_ptr<const DistDvec> ccvec, const double* modata) const override;

    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_rdm12_monomer(std::shared_ptr<const DistDvec> civec, const int i) const override;

    std::shared_ptr<DistDvec> contract_I(std::shared_ptr<const DistDvec> A, std::shared_ptr<Matrix> coef, int offset, int nstA, int nstB, int nstates) const override;
    std::shared_ptr<DistDvec> contract_J(std::shared_ptr<const DistDvec> B, std::shared_ptr<Matrix> coef, int offset, int nstA, int nstB, int nstates) const override;
};

}

#endif
