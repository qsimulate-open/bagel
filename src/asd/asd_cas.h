//
// BAGEL - Parallel electron correlation program.
// Filename: asd_cas.h
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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
