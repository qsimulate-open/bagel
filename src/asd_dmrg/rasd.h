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

namespace bagel {

/// Implementation of ASD_DMRG_Base with RAS model
class RASD : public ASD_DMRG {
  protected:
    std::shared_ptr<DMRG_Block> compute_first_block(std::vector<std::shared_ptr<PTree>> input, std::shared_ptr<const Reference> ref) override;
    std::shared_ptr<DMRG_Block> grow_block(std::vector<std::shared_ptr<PTree>> input, std::shared_ptr<const Reference> ref, std::shared_ptr<DMRG_Block> left, const int site) override;
    std::shared_ptr<DMRG_Block> decimate_block(std::shared_ptr<PTree> input, std::shared_ptr<const Reference> ref, std::shared_ptr<DMRG_Block> system, std::shared_ptr<DMRG_Block> environment, const int site) override;

  public:
    RASD(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer);

  private:
    void read_restricted(std::shared_ptr<PTree> input, const int site) const;

    std::vector<std::shared_ptr<const RASDvec>> diagonalize_site_RDM(const std::vector<std::shared_ptr<ProductRASCivec>>& civecs) const ;
};

}

#endif
