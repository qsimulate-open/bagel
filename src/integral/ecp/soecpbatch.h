//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: soecpbatch.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_INTEGRAL_ECP_SOECPBATCH_H
#define __SRC_INTEGRAL_ECP_SOECPBATCH_H

#include <src/util/parallel/resources.h>
#include <src/molecule/molecule.h>
#include <src/integral/integral.h>
#include <src/integral/ecp/soangularbatch.h>

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
    int asize_final_, asize_;
    bool swap01_;
    size_t size_block_, size_alloc_;
    double* stack_save_;

    bool allocated_here_;
    std::shared_ptr<StackMem> stack_;

    void common_init();
    void get_data(const double* intermediate, double* data) const;

  public:
    SOECPBatch(const std::array<std::shared_ptr<const Shell>,2>& info, const std::shared_ptr<const Molecule> mol,
                   std::shared_ptr<StackMem> = nullptr);
    ~SOECPBatch();

    const double* data() const { return data_; }
    virtual double* data(const int i) override { assert(i == 0); return data_; }

    const double* data1() const { return data1_; }
    const double* data2() const { return data2_; }

    bool swap01() const { return swap01_; }

    void compute() override;

};

}

#endif

