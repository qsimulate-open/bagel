//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rasd.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifndef __ASD_DMRG_RASD_H
#define __ASD_DMRG_RASD_H

#include <src/asd/dmrg/asd_dmrg.h>
#include <src/asd/dmrg/product_civec.h>
#include <src/ci/fci/mofile.h>

namespace bagel {

/// Implementation of ASD_DMRG_Base with RAS model
class RASD : public ASD_DMRG {
  protected:
    std::shared_ptr<DMRG_Block1> compute_first_block(std::vector<std::shared_ptr<PTree>> input, std::shared_ptr<const Reference> ref) override;
    std::shared_ptr<DMRG_Block1> grow_block(std::vector<std::shared_ptr<PTree>> input, std::shared_ptr<const Reference> ref, std::shared_ptr<DMRG_Block1> left, const int site) override;
    std::shared_ptr<DMRG_Block1> decimate_block(std::shared_ptr<PTree> input, std::shared_ptr<const Reference> ref, std::shared_ptr<DMRG_Block1> system, std::shared_ptr<DMRG_Block1> environment, const int site) override;

  public:
    RASD(const std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref);

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
