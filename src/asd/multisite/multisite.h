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

/// Contains references for isolated sites in the ASD_DMRG algorithm.
class MultiSite {
  protected:
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Reference> hf_ref_;
    std::shared_ptr<const Reference> sref_;

    const int nsites_;
    int charge_;
    int nspin_;
    std::vector<int> active_electrons_;
    std::vector<int> active_sizes_;
    std::vector<int> region_sizes_;

  public:
    // Constructors
    MultiSite(std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref, const int nsites);

    int nsites() const { return nsites_; }
    int charge() const { return charge_; }
    int nspin() const { return nspin_; }
    std::vector<int> active_electrons() const { return active_electrons_; }
    std::vector<int> active_sizes() const { return active_sizes_; }

    std::shared_ptr<const Reference> sref() const { return sref_; }

    void compute();

    // Utility functions
    void localize(std::shared_ptr<const PTree> ldata, std::shared_ptr<const Matrix> fock);
    void project_active(std::shared_ptr<const PTree> pdata);
    void set_active();
    void canonicalize(std::shared_ptr<const Matrix> fock);

    std::shared_ptr<Reference> build_reference(const int site, const std::vector<bool> meanfield) const;

    void run_fci() const;
};

}

#endif
