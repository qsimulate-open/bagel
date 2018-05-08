//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_ras.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __ASD_ASD_RAS_H
#define __ASD_ASD_RAS_H

#include <src/asd/asd.h>

namespace bagel {

class ASD_RAS : public ASD<RASDvec> {
  public:
    ASD_RAS(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DimerRAS> cispace, bool rdm = false)
      : ASD<RASDvec>(input, dimer, cispace, rdm) {};

  private:
    std::shared_ptr<RASDvec> form_sigma(std::shared_ptr<const RASDvec> ccvec, std::shared_ptr<const MOFile> jop) const override;
    std::shared_ptr<RASDvec> form_sigma_1e(std::shared_ptr<const RASDvec> ccvec, const double* modata) const override;

    void sigma_aa(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASCivec> sigma, const double* h1, const double* h2) const;
    void sigma_bb(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASCivec> sigma, const double* h1, const double* h2) const;
    void sigma_ab(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASCivec> sigma, const double* h1, const double* h2) const;

    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_rdm12_monomer(std::shared_ptr<const RASDvec> civec, const int i) const override;

    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> compute_rdm12_last_step(std::shared_ptr<const RASCivec> cbra, std::shared_ptr<const RASDvec> dbra) const;
    void sigma_2a(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASDvec> d) const;
    void sigma_2a1(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASDvec> d) const;
    void sigma_2a2(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASDvec> d) const;
    void sigma_2a1_aa(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASDvec> e) const;
    void sigma_2a2_bb(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASDvec> e) const;
    void sigma_2a3_ba(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASDvec> e) const;
    void sigma_2a4_ab(std::shared_ptr<const RASCivec> cc, std::shared_ptr<RASDvec> e) const;

    std::shared_ptr<RASDvec> contract_I(std::shared_ptr<const RASDvec> A, std::shared_ptr<Matrix> coef, int offset, int nstA, int nstB, int nstates) const override;
    std::shared_ptr<RASDvec> contract_J(std::shared_ptr<const RASDvec> B, std::shared_ptr<Matrix> coef, int offset, int nstA, int nstB, int nstates) const override;
};

}

#endif
