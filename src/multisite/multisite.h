//
// BAGEL - Parallel electron correlation program.
// Filename: multisite.h
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

#ifndef __BAGEL_MULTISITE_MULTISITE_H
#define __BAGEL_MULTISITE_MULTISITE_H

#include <memory>
#include <vector>

#include <src/input/input.h>
#include <src/wfn/reference.h>

namespace bagel {

/// Contains references for isolated sites in the ASD_DMRG algorithm.
class MultiSite {
  protected:
    std::shared_ptr<const PTree> input_;

    std::vector<std::shared_ptr<const Reference>> isolated_refs_; ///< Reference objects of the isolated monomers BEFORE active spaces have been chosen
    std::vector<std::shared_ptr<const Reference>> active_refs_;   ///< Reference objects of the isolated monomers AFTER the active spaces have been chosen

    std::vector<std::shared_ptr<const Geometry>> geoms_; ///< hold onto original geometry objects for nbasis and natom information

    std::shared_ptr<const Reference> sref_; ///< Super-reference, i.e., Reference of whole multisite system
    double active_thresh_;                  ///< overlap threshold for inclusion in the active space

    std::vector<std::pair<int, int>> closed_bounds_; ///< list of [start, end) pairs for the closed spaces of each site
    std::vector<std::pair<int, int>> active_bounds_; ///< list of [start, end) pairs for the active spaces of each site
    std::vector<std::pair<int, int>> virt_bounds_; ///< list of [start, end) pairs for the virtual spaces of each site
    std::vector<std::pair<int, int>> occ_act_bounds_; ///< list of [start, end) pairs for all of the occupied active orbitals (orbitals occupied in a HF sense)

    const int nsites_;

  public:
    // Constructors
    MultiSite(std::shared_ptr<const PTree> input, std::vector<std::shared_ptr<const Reference>> refs); ///< Conjoins the provided Reference objects

    // Return functions
    std::vector<std::shared_ptr<const Reference>> isolated_refs() const { return isolated_refs_; }
    std::vector<std::shared_ptr<const Reference>> active_refs() const { return active_refs_; }

    std::shared_ptr<const Reference> isolated_refs(const int i) const { return isolated_refs_.at(i); }
    std::shared_ptr<const Reference> active_refs(const int i) const { return active_refs_.at(i); }

    std::shared_ptr<const Reference> conv_to_ref() const { return sref_; }

    // Utility functions
    /// Sets active space of sref_ using overlaps with isolated_ref_ active spaces
    void set_active(std::shared_ptr<const PTree> idata);
    /// Localizes active space and uses the given Fock matrix to diagonalize the subspaces
    void localize(std::shared_ptr<const PTree> idata, std::shared_ptr<const Matrix> fock);

    void scf(std::shared_ptr<const PTree> idata); ///< Driver for preparation of sites for ASD_DMRG

    /// Creates a Reference object for an ASD calculation
    std::shared_ptr<Reference> build_reference(const int site, const std::vector<bool> meanfield) const;
};

}

#endif
