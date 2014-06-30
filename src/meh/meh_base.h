//
// BAGEL - Parallel electron correlation program.
// Filename: meh_base.h
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

#ifndef __MEH_MEH_BASE_H
#define __MEH_MEH_BASE_H

#include <src/dimer/dimer.h>
#include <src/dimer/dimer_jop.h>
#include <src/dimer/dimer_prop.h>
#include <src/meh/meh_spin.h>

namespace bagel {

/// Enumeration of the possible couplings between dimer blocks
enum class Coupling {
  none = 0,       ///< no coupling
  diagonal = 1,   ///< no change in occupation patterns
  aET = 2,        ///< alpha transfer (A --> B)
  inv_aET = -2,   ///< inverse alpha transfer (B --> A)
  bET = 3,        ///< beta transfer (A --> B)
  inv_bET = -3,   ///< inverse beta transfer (B --> A)
  abFlip = 4,     ///< alpha-beta flip
  baFlip = -4,    ///< beta-alpha flip
  abET = 5,       ///< alpha+beta transfer (A --> B)
  inv_abET = -5,  ///< inverse alpha+beta (B --> A)
  aaET = 6,       ///< alpha+alpha transfer (A --> B)
  inv_aaET = -6,  ///< inverse alpha+alpha transfer (B --> A)
  bbET = 7,       ///< beta+beta transfer (A --> B)
  inv_bbET = -7   ///< inverse beta+beta transfer (B --> A)
};

/// Specifies a single block of a model Hamiltonian
struct ModelBlock {
  std::pair<int, int> S_;
  std::pair<int, int> charge_;
  std::pair<int, int> M_;
  int nstates_;

  ModelBlock(std::pair<int, int> S, std::pair<int, int> q, std::pair<int, int> M, const int nstates) :
    S_(S), charge_(q), M_(M), nstates_(nstates) {}
};

/// Base class for MEH (to be renamed ASD). Contains everything that doesn't need to be templated.
class MEH_base {
  protected:
    std::shared_ptr<const Dimer> dimer_;

    std::shared_ptr<DimerJop> jop_;

    std::shared_ptr<Matrix> hamiltonian_; ///< if stored_ is true, Hamiltonian is stored here
    std::shared_ptr<Matrix> adiabats_; ///< Eigenvectors of adiabatic states
    std::vector<std::pair<std::string, std::shared_ptr<Matrix>>> properties_;

    std::unique_ptr<double[]> denom_; ///< denominator used for Davidson

    int max_spin_;
    std::shared_ptr<MEHSpin> spin_; ///< Sparsely stored spin matrix

    std::vector<double> energies_; ///< Adiabatic energies

    // Total system quantities
    int dimerstates_; ///< Total size of dimer Hamiltonian. Counted up during initialization

    // Options
    int nstates_;
    int nspin_;
    int charge_;
    int max_iter_;
    int nguess_;
    int davidson_subspace_;

    bool store_matrix_;
    bool dipoles_;

    double thresh_;
    double print_thresh_;

    std::vector<std::vector<ModelBlock>> models_to_form_; ///< Contains specifications to construct model spaces
    std::vector<std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>>> models_; ///< models that have been built

  public:
    MEH_base(const std::shared_ptr<const PTree> input, std::shared_ptr<const Dimer> dimer);

    virtual void compute() = 0;
    virtual void compute_rdm() const = 0;

    std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> model(const int i) { return models_[i]; }

    std::vector<double> energy() const { return energies_; }
    double energy(const int i) const { return energies_.at(i); }

};

}

#endif
