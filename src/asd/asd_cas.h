//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_cas.h
// Copyright (C) 2012 Toru Shiozaki
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

#ifndef __ASD_ASD_CAS_H
#define __ASD_ASD_CAS_H

#include <src/asd/asd.h>

namespace bagel {

class ASD_CAS : public ASD<CASDvec> {
  public:
    ASD_CAS(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DimerCAS> cispace, bool rdm = false)
      : ASD<CASDvec>(input, dimer, cispace, rdm) {
      cispace_->intermediates();
    }

  private:
    std::shared_ptr<CASDvec> form_sigma(std::shared_ptr<const CASDvec> ccvec, std::shared_ptr<const MOFile> jop) const override;
    std::shared_ptr<CASDvec> form_sigma_1e(std::shared_ptr<const CASDvec> ccvec, const double* modata) const override;

    void sigma_aa(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, const double* h1, const double* h2) const;
    void sigma_2ab_1(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d) const;
    void sigma_2ab_2(std::shared_ptr<const Dvec> d, std::shared_ptr<Dvec> e, const double* mo2e_ptr) const;
    void sigma_2ab_3(std::shared_ptr<Civec> sigma, std::shared_ptr<Dvec> e) const;

    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_rdm12_monomer(std::shared_ptr<const CASDvec> civec, const int i) const override;

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>> compute_rdm12_last_step(std::shared_ptr<const Civec>, std::shared_ptr<const Dvec>) const;
    void sigma_2a1(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d) const;
    void sigma_2a2(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d) const;

    std::shared_ptr<CASDvec> contract_I(std::shared_ptr<const CASDvec> A, std::shared_ptr<Matrix> coef, int offset, int nstA, int nstB, int nstates) const override;
    std::shared_ptr<CASDvec> contract_J(std::shared_ptr<const CASDvec> B, std::shared_ptr<Matrix> coef, int offset, int nstA, int nstB, int nstates) const override;
};

}

#endif
