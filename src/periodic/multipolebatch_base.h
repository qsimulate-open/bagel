//
// BAGEL - Parallel electron correlation program.
// Filename: multipolebatch_base.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_PERIODIC_MULTIPOLEBATCH_BASE_H
#define __SRC_PERIODIC_MULTIPOLEBATCH_BASE_H

#include <src/parallel/resources.h>
#include <src/molecule/atom.h>
#include <src/integral/integral.h>
#include <src/periodic/multipole.h>

namespace bagel {

class MultipoleBatch_base : public Integral {
  protected:
    std::array<std::shared_ptr<const Shell>,2> basisinfo_;
    std::shared_ptr<const Atom> site_;
    bool spherical_;

    bool swap01_;
    double* P_;
    std::array<double, 3> AB_;

    size_t size_block_, size_alloc_, size_allocated_;
    int primsize_, prim0size_, prim1size_;
    int contsize_, cont0size_, cont1size_;
    int asize_, size_final_, amax_, amax1_;
    int num_multipoles_;
    std::complex<double>* multipole_;

    bool allocated_here_;
    std::shared_ptr<StackMem> stack_;

    double* data_;
    double* buff_, stack_save_;
    std::complex<double>* multipole_buff;

    int amapping_[ANG_VRR_END * ANG_VRR_END * ANG_VRR_END];

    void init();
    void allocate_arrays(const size_t ps);
    void compute_ss(const double thresh);
    void allocate_data(const int asize_final_sph);

  public:
    MultipoleBatch_base(const std::array<std::shared_ptr<const Shell>,2>& shells, const std::shared_ptr<const Atom> atom,
                        std::shared_ptr<StackMem> stack = nullptr);
    ~MultipoleBatch_base() { }

    std::shared_ptr<const Atom> site() const { return site_; }

    virtual void compute() = 0;

    double* data(const int i) override { assert(i == 0); return data_; }
    const double* data() const { return data_; }

    bool swap01() const { return swap01_; }

};

}

#endif
