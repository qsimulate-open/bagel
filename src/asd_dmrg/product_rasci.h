//
// BAGEL - Parallel electron correlation program.
// Filename: asd_dmrg/product_rasci.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __BAGEL_ASD_DMRG_PRODUCT_RASCI_H
#define __BAGEL_ASD_DMRG_PRODUCT_RASCI_H

#include <src/dimer/dimer_jop.h>

#include <src/asd_dmrg/product_civec.h>
#include <src/asd_dmrg/dmrg_block.h>
#include <src/asd_dmrg/block_operators.h>

namespace bagel {

/// Diagonalizes a CI space \f$|L\rangle\otimes|\mbox{RAS}\rangle\f$ where \f$|L\rangle\f$ is a set of (DMRG) block states and \f$|\mbox{RAS}\rangle\f$ is a normal RAS wavefunction
class ProductRASCI {
  protected:
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Reference> ref_;
    std::shared_ptr<const DMRG_Block> left_;
    std::shared_ptr<const BlockOperators> blockops_;

    int max_iter_; ///< maximum number of iterations to perform
    int davidson_subspace_; ///< maximum number of trial vectors to keep during the Davidson iterations
    int nguess_; ///< number of guess vectors to start

    double thresh_; ///< convergence threshold
    double print_thresh_; ///< printing threshold

    int batchsize_; ///< batchsize used in sigma_aa and sigma_bb portions of FormSigma

    int nelea_; ///< total number of \f$\alpha\f$ electrons in product space
    int neleb_; ///< total number of \f$\alpha\f$ electrons in product space

    int ncore_; ///< number of core/closed orbitals on single site
    int norb_;  ///< number of active orbitals on single site, \f$\mbox{norb}=\sum_i\mbox{ras[i]}\f$
    std::array<int,3> ras_; ///< {RASI, RASII, RASIII} orbitals on single site
    int max_holes_; ///< maximum number of holes in RASI on single site
    int max_particles_; ///< maximum number of particles in RASIII on single site

    int nstates_; ///< number of roots to compute for the product space
    std::vector<double> energy_; ///< total energies of the roots computed

    std::vector<std::shared_ptr<ProductRASCivec>> cc_; ///< CI vector at convergence

    std::shared_ptr<DimerJop> jop_; ///< MO integrals
    std::shared_ptr<RASSpace> space_; ///< space of determinants used in calculation
    std::shared_ptr<ProductRASCivec> denom_; ///< denominator used in Davidson routine

  public:
    /// Constructor
    ProductRASCI(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref, std::shared_ptr<const DMRG_Block> left);

    void compute();

    // returns members
    int norb() const { return norb_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
    int ncore() const { return ncore_; }
    int ras(const int i) const { return ras_[i]; }
    double core_energy() const { return jop_->core_energy(); }

    std::shared_ptr<const DimerJop> jop() const { return jop_; }

    std::vector<double> energy() const { return energy_; }
    double energy(const int i) const { return energy_.at(i); }

    const std::vector<std::shared_ptr<ProductRASCivec>>& civectors() const { return cc_; }

  private:
    void construct_denom(); ///< construct denominator

    void common_init(); ///< initializer
    void model_guess(std::vector<std::shared_ptr<ProductRASCivec>>& out);
    std::vector<std::pair<std::bitset<nbit__>, std::bitset<nbit__>>> detseeds(const int ndet); ///< generate spin-adapted guess configurations

    void print_header() const; ///< print functions
};

}

#endif

