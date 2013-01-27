//
// BAGEL - Parallel electron correlation program.
// Filename: osint.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __src_osint_osint_h
#define __src_osint_osint_h

#include <stddef.h>
#include <array>
#include <vector>
#include <src/scf/shell.h>
#include <src/rysint/sortlist.h>
#include <array>
#include <memory>
#include <src/parallel/resources.h>
#include <src/rysint/integral.h>

namespace bagel {

class OSInt : public Integral {
  protected:

    std::array<std::shared_ptr<const Shell>,2> basisinfo_;
    bool spherical_;

    double* data_;
    std::vector<double> xp_, xa_, xb_, rho_, p_;
    std::vector<double> coeffsx_, coeffsy_, coeffsz_;
    std::vector<double> coefftx_, coeffty_, coefftz_;
    std::array<double,3> AB_;

    int ang0_, ang1_, cont0_, cont1_, prim0_, prim1_;

    int amax_, amax1_, amin_, asize_, asize_final_, asize_intermediate_;
    std::vector<int> amapping_;

    bool swap01_;

    const SortList sort_;

    const int deriv_rank_;
    size_t size_alloc_;
    size_t size_block_;

    double* stack_save_;

    virtual void perform_VRR(double*) {};
    void perform_contraction(const int, const double*, const int, const int, double*,
                             const std::vector<std::vector<double> >&, const std::vector<std::pair<int, int> >&, const int,
                             const std::vector<std::vector<double> >&, const std::vector<std::pair<int, int> >&, const int);

    bool allocated_here_;
    std::shared_ptr<StackMem> stack_;

  public:
    // deriv rank negative means multipole integrals
    OSInt(const std::array<std::shared_ptr<const Shell>,2>&, const int deriv = 0, std::shared_ptr<StackMem> = std::shared_ptr<StackMem>());
    ~OSInt();

    constexpr static int nblocks() { return 1; }
    double* data(const int i) override { return data_; }
    const double* data() const { return data_; };

    // since this is convenient for gradient evaluation....
    bool swap01() const { return swap01_; };

    size_t size_block() const { return size_block_; };
};

}

#endif

