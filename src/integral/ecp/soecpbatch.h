//
// BAGEL - Parallel electron correlation program.
// Filename: soecpbatch.h
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


#ifndef __SRC_INTEGRAL_ECP_SOECPBATCH_H
#define __SRC_INTEGRAL_ECP_SOECPBATCH_H

#include <src/parallel/resources.h>
#include <src/molecule/molecule.h>
#include <src/integral/integral.h>
#include <src/integral/ecp/sobatch.h>

namespace bagel {

class SOECPBatch : public Integral {
  protected:

    int max_iter_;
    double integral_thresh_;

    std::array<std::shared_ptr<const Shell>,2> basisinfo_;
    std::shared_ptr<const Molecule> mol_;

    bool spherical_;

    double* data_;
    double* data1_;
    double* data2_;

    int ang0_, ang1_, cont0_, cont1_;
    int amax_, amax1_, amin_, asize_, asize_final_, asize_intermediate_;
    bool swap01_;
    size_t size_alloc_;
    double* stack_save_;

    bool allocated_here_;
    std::shared_ptr<StackMem> stack_;

    void common_init();
    void get_data(double* intermediate, double* data);

  public:
    SOECPBatch(const std::array<std::shared_ptr<const Shell>,2>& info, const std::shared_ptr<const Molecule> mol,
                   std::shared_ptr<StackMem> = nullptr);
    ~SOECPBatch();

    double* data() { return data_; }
    virtual double* data(const int i) override { assert(i == 0); return data_; }

    double* data1() { return data1_; }
    double* data2() { return data2_; }

    bool swap01() const { return swap01_; }

    void compute() override;

};

}

#endif

