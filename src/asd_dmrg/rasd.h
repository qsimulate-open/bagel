//
// BAGEL - Parallel electron correlation program.
// Filename: rasd.h
// Copyright (C) 2014 Shane Parker
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

#ifndef __ASD_DMRG_RASD_H
#define __ASD_DMRG_RASD_H

#include <src/asd_dmrg/asd_dmrg.h>
#include <src/asd_dmrg/product_civec.h>
#include <src/fci/mofile.h>

namespace bagel {

/// Implementation of ASD_DMRG_Base with RAS model
class RASD : public ASD_DMRG {
  protected:
    std::shared_ptr<DMRG_Block1> compute_first_block(std::vector<std::shared_ptr<PTree>> input, std::shared_ptr<const Reference> ref) override;
    std::shared_ptr<DMRG_Block1> grow_block(std::vector<std::shared_ptr<PTree>> input, std::shared_ptr<const Reference> ref, std::shared_ptr<DMRG_Block1> left, const int site) override;
    std::shared_ptr<DMRG_Block1> decimate_block(std::shared_ptr<PTree> input, std::shared_ptr<const Reference> ref, std::shared_ptr<DMRG_Block1> system, std::shared_ptr<DMRG_Block1> environment, const int site) override;

  public:
    RASD(const std::shared_ptr<const PTree> input, std::shared_ptr<MultiSite> multisite);

  private:
    void read_restricted(std::shared_ptr<PTree> input, const int site) const;

    /// diagonalize site RDM
    std::map<BlockKey, std::shared_ptr<const RASDvec>> diagonalize_site_RDM(const std::vector<std::shared_ptr<ProductRASCivec>>& civecs, const double perturbation = 0.0) const;
    /// diagonalize site+block RDM
    std::map<BlockKey, std::vector<std::shared_ptr<ProductRASCivec>>> diagonalize_site_and_block_RDM(const std::vector<std::shared_ptr<ProductRASCivec>>& civecs, const double perturbation = 0.0) const;

    /// compute 2e hamiltonian for RASCivecs
    std::shared_ptr<Matrix> compute_sigma2e(std::shared_ptr<const RASDvec> cc, std::shared_ptr<const MOFile> jop) const;
    /// compute 2e hamiltonian for ProductRASCivecs
    std::shared_ptr<Matrix> compute_sigma2e(const std::vector<std::shared_ptr<ProductRASCivec>>& cc, std::shared_ptr<const DimerJop> jop) const;
    // compute spin matrix for RASCivecs
    std::shared_ptr<Matrix> compute_spin(std::shared_ptr<const RASDvec> cc) const;
    // compute spin matrix for ProductRASCivecs
    std::shared_ptr<Matrix> compute_spin(std::vector<std::shared_ptr<ProductRASCivec>> cc) const;

    void apply_perturbation(std::shared_ptr<const RASBlockVectors> cc, std::vector<GammaSQ> oplist, std::map<BlockKey, std::shared_ptr<const RASDeterminants>>& detmap,
                            const double weight, std::map<BlockKey, std::vector<std::tuple<double, std::shared_ptr<const Matrix>>>>& outer_products) const;
    void apply_perturbation(const std::vector<std::shared_ptr<ProductRASCivec>>& ccvec, const BlockKey cckey, std::vector<GammaSQ> oplist,
                        const double weight, std::map<BlockKey, std::vector<std::tuple<double, std::vector<std::shared_ptr<ProductRASCivec>>>>>& outer_products) const;
};

}

#endif
