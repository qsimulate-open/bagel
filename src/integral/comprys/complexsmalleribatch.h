//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexsmalleribatch.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_INTEGRAL_COMPRYS_COMPLEXSMALLERIBATCH_H
#define __SRC_INTEGRAL_COMPRYS_COMPLEXSMALLERIBATCH_H

#include <src/integral/comprys/complexeribatch.h>

namespace bagel {

class ComplexSmallERIBatch : public Integral_base<std::complex<double>> {
  protected:
    std::complex<double>* data_;

    // size info
    size_t size_block_;
    size_t size_alloc_;

    const std::array<std::shared_ptr<const Shell>,3> shells_;

    std::shared_ptr<StackMem> stack_;

  public:
    ComplexSmallERIBatch(std::array<std::shared_ptr<const Shell>,4> info, const double dummy);

    ~ComplexSmallERIBatch();

    void compute();

    std::complex<double>* data(const int i) override { return data_+i*size_block_; }

    void eri_compute(std::complex<double>* eri) const;

    size_t size_block() const { return size_block_; }
    size_t size_alloc() const { return size_alloc_; }

    // need to return 6 blocks
    constexpr static int Nblocks() { return 6; }
};

}

#endif
