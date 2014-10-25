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
#include <src/ras/sparse_ij.h>
#include <src/asd_dmrg/phi_k_lists.h>
#include <src/asd_dmrg/phi_ijk_lists.h>

namespace bagel {

class FormSigmaProdRAS {
  protected:
    int batchsize_; ///< batchsize used in \f$\alpha\alpha\f$ and \f$\beta\beta\f$ parts of pure RAS

  public:
    FormSigmaProdRAS(const int b = 512) : batchsize_(b) {}

    /// Applies Hamiltonian to cc using the provided MOFile, skipping the vectors marked as converged
    std::vector<std::shared_ptr<ProductRASCivec>> operator()(const std::vector<std::shared_ptr<ProductRASCivec>>& ccvec, std::shared_ptr<const BlockOperators> blockops, std::shared_ptr<DimerJop> jop, const std::vector<bool>& conv) const;
    std::vector<std::shared_ptr<ProductRASCivec>> diagonal(const std::vector<std::shared_ptr<ProductRASCivec>>& ccvec, std::shared_ptr<const BlockOperators> blockops, std::shared_ptr<DimerJop> jop, const std::vector<bool>& conv) const;

  private:
    // Helper functions for sigma formation
    void pure_block_and_ras(std::shared_ptr<const ProductRASCivec> cc, std::shared_ptr<ProductRASCivec> sigma, std::shared_ptr<const BlockOperators> blockops, std::shared_ptr<DimerJop> jop) const;
    void interaction_terms(std::shared_ptr<const ProductRASCivec> cc, std::shared_ptr<ProductRASCivec> sigma, std::shared_ptr<const BlockOperators> blockops, std::shared_ptr<DimerJop> jop) const;
    void diagonal_terms(std::shared_ptr<const ProductRASCivec> cc, std::shared_ptr<ProductRASCivec> sigma, std::shared_ptr<const BlockOperators> blockops, std::shared_ptr<DimerJop> jop) const;

    /// Branch 1: \f$\alpha^\dagger, \alpha^\dagger\alpha^\dagger\f$
    void aET_branch(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 2: \f$\beta^\dagger, \alpha^\dagger\beta^\dagger\f$
    void bET_branch(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 3: \f$\alpha, \alpha\alpha \f$
    void aHT_branch(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 4: \f$\beta, \alpha\beta, \beta\beta\f$
    void bHT_branch(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 5: \f$\alpha^\dagger\alpha\f$
    void aexc_branch(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 6: \f$\beta^\dagger\beta\f$
    void bexc_branch(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 7: \f$\alpha^\dagger\beta\f$
    void abflip_branch(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;
    /// Branch 8: \f$\beta^\dagger\alpha\f$
    void baflip_branch(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blocksops) const;

    /// Computes 3-operator aET terms
    void compute_sigma_3aET(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blockops,  std::shared_ptr<DimerJop> jop) const;

    /// Computes 3-operator aHT terms
    void compute_sigma_3aHT(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blockops,  std::shared_ptr<DimerJop> jop) const;

    /// Computes 3-operator bET terms
    void compute_sigma_3bET(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blockops,  std::shared_ptr<DimerJop> jop) const;

    /// Computes 3-operator bHT terms
    void compute_sigma_3bHT(std::shared_ptr<const RASBlockVectors> cc, std::shared_ptr<ProductRASCivec> sigma_sector, std::shared_ptr<const BlockOperators> blockops,  std::shared_ptr<DimerJop> jop) const;

    /// Computes \f$\hat S_p = \sum_{ijk} i^\dagger_\alpha j^\dagger_\alpha k_\alpha (jk|pi)\f$ or \f$\hat S_p^\dagger\f$,
    /// depending on what type of PhiIJKLists is given to it.
    void resolve_S_aaa(const RASBlockVectors& cc, RASBlockVectors& sigma, const double* Jp, const PhiIJKLists& phi_ijk) const;

    /// Computes \f$\hat S_p = \sum_{i,j,k} i^\dagger_\alpha j^\dagger_\alpha k_\alpha (jk|pi)f$
    void resolve_S_adag_adag_a(const RASCivecView cc, RASCivecView sigma, const double* Jp, std::shared_ptr<const PhiIJKLists> phi_ijk = nullptr) const;

    /// Computes \f$\hat S_p = \sum_{i,j,k} i^\dagger_\alpha j_\alpha k_\alpha (ij|pk)f$
    void resolve_S_adag_a_a(const RASCivecView cc, RASCivecView sigma, const double* Jp, std::shared_ptr<const PhiIJKLists> phi_ijk = nullptr) const;

    /// Computes \f$\hat S_p = \sum_{i,j,k} i^\dagger_\beta j_\beta k_\alpha (ij|pk)f$ and \f$\hat S_p = \sum_{i,j,k}
    void resolve_S_abb(const RASCivecView cc, RASCivecView sigma, const double* Jp, std::shared_ptr<PhiKLists> phik = nullptr, std::shared_ptr<Sparse_IJ> sparse = nullptr) const;
};

}

#endif

