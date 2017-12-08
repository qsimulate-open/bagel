//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multisite.h
// Copyright (C) 2014 Shane Parker
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

#ifndef __BAGEL_MULTISITE_MULTISITE_H
#define __BAGEL_MULTISITE_MULTISITE_H

#include <memory>
#include <vector>
#include <src/util/input/input.h>
#include <src/wfn/reference.h>

namespace bagel {

// MultiSite provides system information for ASD-DMRG solver
class MultiSite {
  protected:
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Reference> hf_ref_;
    std::shared_ptr<const Reference> sref_;

    // system info
    const int nsites_;
    int charge_;
    int nspin_;
    std::vector<int> active_electrons_; // specify number of electrons on each site as initial guess
    std::vector<int> active_sizes_;     // number of active orbitals localized on each site
    std::vector<int> region_sizes_;     // number of atoms on each site; atoms should be ordered by sites in molecular geometry

  public:
    // constructor
    MultiSite(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref, const int nsites);

    void compute();

    // utility functions
    void localize(std::shared_ptr<const PTree> ldata, std::shared_ptr<const Matrix> fock);
    void set_active_orbitals();
    void canonicalize(std::shared_ptr<const Matrix> fock);

    // return functions
    int nsites() const { return nsites_; }
    int charge() const { return charge_; }
    int nspin() const { return nspin_; }
    std::vector<int> active_electrons() const { return active_electrons_; }
    std::vector<int> active_sizes() const { return active_sizes_; }
    std::shared_ptr<const Reference> conv_to_ref() const { return sref_; }
    
    // prepare input reference for ASD-DMRG
    std::shared_ptr<Reference> build_reference(const int site, const std::vector<bool> meanfield) const;
};

}

#endif
