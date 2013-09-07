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


#ifndef __SRC_INTEGRAL_OS_OSINT_H
#define __SRC_INTEGRAL_OS_OSINT_H

#include <src/molecule/shell.h>
#include <src/integral/sortlist.h>
#include <src/math/xyzfile.h>
#include <src/parallel/resources.h>
#include <src/integral/integral.h>

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

    size_t size_alloc_;
    size_t size_block_;

    double* stack_save_;

    virtual void perform_VRR(double*) { throw std::logic_error("OSInt::perform_VRR should not be called"); }
    void perform_contraction(const int, const double*, const int, const int, double*,
                             const std::vector<std::vector<double>>&, const std::vector<std::pair<int, int>>&, const int,
                             const std::vector<std::vector<double>>&, const std::vector<std::pair<int, int>>&, const int);

    bool allocated_here_;
    std::shared_ptr<StackMem> stack_;

    virtual int nblocks() const = 0;
    virtual int nrank() const = 0;
    void common_init();

  public:
    // deriv rank negative means multipole integrals
    OSInt(const std::array<std::shared_ptr<const Shell>,2>&, std::shared_ptr<StackMem> = std::shared_ptr<StackMem>());
    ~OSInt();

    double* data(const int i) override { assert(i < nblocks()); return data_+i*size_block_; }
    const double* data() const { return data_; }

    // since this is convenient for gradient evaluation....
    bool swap01() const { return swap01_; }

    size_t size_block() const { return size_block_; }

    virtual std::shared_ptr<GradFile> compute_gradient(std::shared_ptr<const Matrix> d, const int iatom0, const int iatom1, const int natom) const;
};

}

#endif

