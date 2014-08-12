//
// BAGEL - Parallel electron correlation program.
// Filename: asd_dmrg/form_sigma.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __BAGEL_ASD_DMRG_FORM_SIGMA_H
#define __BAGEL_ASD_DMRG_FORM_SIGMA_H

#include <src/dimer/dimer_jop.h>
#include <src/asd_dmrg/product_civec.h>
#include <src/asd_dmrg/block_operators.h>

namespace bagel {

class FormSigmaProdRAS {
  protected:
    int batchsize_; ///< batchsize used in \f$\alpha\alpha\f$ and \f$\beta\beta\f$ parts of pure RAS

  public:
    FormSigmaProdRAS(const int b = 512) : batchsize_(b) {}

    /// Applies Hamiltonian to cc using the provided MOFile, skipping the vectors marked as converged
    std::vector<std::shared_ptr<ProductRASCivec>> operator()(std::vector<std::shared_ptr<ProductRASCivec>>& ccvec, std::shared_ptr<const BlockOperators> blockops, std::shared_ptr<const DimerJop> jop, const std::vector<bool>& conv) const;

  private:
    // Helper functions for sigma formation
    void pure_block_and_ras(std::shared_ptr<const ProductRASCivec> cc, std::shared_ptr<ProductRASCivec> sigma, std::shared_ptr<const BlockOperators> blockops, std::shared_ptr<const DimerJop> jop) const;
    void interaction_terms(std::shared_ptr<const ProductRASCivec> cc, std::shared_ptr<ProductRASCivec> sigma, std::shared_ptr<const BlockOperators> blockops, std::shared_ptr<const DimerJop> jop) const;

    /// Branch 1: \f$\alpha^\dagger, \alpha^\dagger\alpha^\dagger\f$
    void branch_1(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 2: \f$\beta^\dagger, \alpha^\dagger\beta^\dagger\f$
    void branch_2(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 3: \f$\alpha, \alpha\alpha, \alpha^\dagger\alpha\alpha, \beta^\dagger \alpha\alpha, (\alpha^\dagger\beta\alpha)\f$
    void branch_3(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 4: \f$\beta, \alpha\beta, \beta\beta, (\beta^\dagger\alpha\beta), \beta^\dagger\beta\beta, (\alpha^\dagger\alpha\beta)\f$
    void branch_4(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 5: \f$\alpha^\dagger\alpha, \alpha^\dagger\alpha^\dagger\alpha, (\beta^\dagger\alpha^\dagger\alpha)\f$
    void branch_5(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 6: \f$\beta^\dagger\beta, \beta^\dagger\beta^\dagger\beta, (\alpha^\dagger\beta^\dagger\beta)\f$
    void branch_6(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 7: \f$\alpha^\dagger\beta, (\beta^\dagger\alpha^\dagger\beta)\f$
    void branch_7(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 8: \f$\beta^\dagger\alpha, (\alpha^\dagger\beta^\dagger\alpha)\f$
    void branch_8(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
};

}

#endif

