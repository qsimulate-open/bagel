//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg/product_rasci.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#ifndef __BAGEL_ASD_DMRG_PRODUCT_RASCI_H
#define __BAGEL_ASD_DMRG_PRODUCT_RASCI_H

#include <src/asd/dimer/dimer_jop.h>
#include <src/asd/dmrg/product_civec.h>
#include <src/asd/dmrg/dmrg_block.h>
#include <src/asd/dmrg/block_operators.h>

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

    bool preconverge_; ///< whether to pre-diagonalize the initial guess using just the "diagonal" part of the hamiltonian
    double preconv_thresh_; ///< threshold at which to preconverge with diagonal portions
    int preconv_iter_; ///< maximum number of iterations for preconvergence

    int batchsize_; ///< batchsize used in sigma_aa and sigma_bb portions of FormSigma

    int nelea_; ///< total number of \f$\alpha\f$ electrons in product space
    int neleb_; ///< total number of \f$\alpha\f$ electrons in product space

    int ncore_; ///< number of core/closed orbitals on single site
    int norb_;  ///< number of active orbitals on single site, \f$\mbox{norb}=\sum_i\mbox{ras[i]}\f$
    std::array<int,3> ras_; ///< {RASI, RASII, RASIII} orbitals on single site
    int max_holes_; ///< maximum number of holes in RASI on single site
    int max_particles_; ///< maximum number of particles in RASIII on single site

    int nstate_; ///< number of roots to compute for the product space
    std::vector<double> energy_; ///< total energies of the roots computed

    std::vector<std::shared_ptr<ProductRASCivec>> cc_; ///< CI vector at convergence

    std::shared_ptr<DimerJop> jop_; ///< MO integrals
    std::shared_ptr<RASSpace> space_; ///< space of determinants used in calculation
    std::shared_ptr<ProductRASCivec> denom_; ///< denominator used in Davidson routine

    double max_coulomb_site_; ///< maximum \f$(ij|ij)^{1/2}\f$ where \f$i\f$ and \f$j\f$ are site orbitals, used for screening
    std::shared_ptr<Matrix> site_block_coulomb_; ///< \f$M(i,p) = (ip|ip)^{1/2}\f$ where \f$i\in\mbox{isite}\f$ and \f$p\in\mbox{block}\f$

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

    std::shared_ptr<DimerJop> jop() const { return jop_; }

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

