//
// BAGEL - Parallel electron correlation program.
// Filename: meh_cas.h
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

#ifndef __MEH_MEH_CAS_H
#define __MEH_MEH_CAS_H

#include <src/meh/meh.h>

namespace bagel {

class MEH_CAS : public MultiExcitonHamiltonian<Dvec> {
   public:
      MEH_CAS(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DimerCAS> cispace) :
        MultiExcitonHamiltonian<Dvec>(input, dimer, cispace)
      {
        cispace_->intermediates();
      }

   private:
      std::shared_ptr<Dvec> form_sigma(std::shared_ptr<const Dvec> ccvec, std::shared_ptr<const MOFile> jop) const override;
      std::shared_ptr<Dvec> form_sigma_1e(std::shared_ptr<const Dvec> ccvec, const double* modata) const override;

      void sigma_aa(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, const double* h1, const double* h2) const;
      void sigma_2ab_1(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d) const;
      void sigma_2ab_2(std::shared_ptr<Dvec> d, std::shared_ptr<Dvec> e, const double* mo2e_ptr) const;
      void sigma_2ab_3(std::shared_ptr<Civec> sigma, std::shared_ptr<Dvec> e) const;
};

}

#endif
