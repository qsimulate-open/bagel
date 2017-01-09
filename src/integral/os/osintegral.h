//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: osintegral.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_INTEGRAL_OS_OSINTEGRAL_H
#define __SRC_INTEGRAL_OS_OSINTEGRAL_H

#include <src/molecule/shell.h>
#include <src/integral/sortlist.h>
#include <src/util/math/xyzfile.h>
#include <src/util/parallel/resources.h>
#include <src/integral/integral.h>

namespace bagel {

template <typename DataType, Int_t IntType = Int_t::Standard>
class OSIntegral : public Integral_base<DataType> {
  protected:

    std::array<std::shared_ptr<const Shell>,2> basisinfo_;
    bool spherical_;

    DataType* data_;
    std::vector<double> xp_, xa_, xb_, rho_;
    std::vector<DataType> P_;
    std::vector<DataType> coeffsx_, coeffsy_, coeffsz_;
    std::vector<DataType> coefftx_, coeffty_, coefftz_;
    std::array<double,3> AB_;

    int ang0_, ang1_, cont0_, cont1_, prim0_, prim1_;

    int amax_, amax1_, amin_, asize_, asize_final_, asize_intermediate_;
    std::vector<int> amapping_;

    bool swap01_;

    size_t size_alloc_;
    size_t size_block_;

    DataType* stack_save_;

    virtual void perform_VRR(DataType*) { throw std::logic_error("OSInt::perform_VRR should not be called"); }
    void perform_contraction(const int, const DataType*, const int, const int, DataType*,
                             const std::vector<std::vector<double>>&, const std::vector<std::pair<int, int>>&, const int,
                             const std::vector<std::vector<double>>&, const std::vector<std::pair<int, int>>&, const int);
    virtual DataType get_P(const double coord1, const double coord2, const double exp1, const double exp2, const double one12, const int dim, const bool swap) {
      return (coord1*exp1 + coord2*exp2) * one12;
    }

    bool allocated_here_;
    std::shared_ptr<StackMem> stack_;

    virtual int nblocks() const = 0;
    virtual int nrank() const = 0;
    void common_init();

  public:
    // deriv rank negative means multipole integrals
    OSIntegral(const std::array<std::shared_ptr<const Shell>,2>&, std::shared_ptr<StackMem> = nullptr);
    ~OSIntegral();

    DataType* data(const int i) override { assert(i < nblocks()); return data_+i*size_block_; }
    const DataType* data() const { return data_; }

    // since this is convenient for gradient evaluation....
    bool swap01() const { return swap01_; }

    size_t size_block() const { return size_block_; }
    size_t asize_final() const { return asize_final_; }

    virtual std::shared_ptr<GradFile> compute_gradient(std::shared_ptr<const Matrix> d, const int iatom0, const int iatom1, const int natom) const;
};

using OSInt = OSIntegral<double,Int_t::Standard>;

}

extern template class bagel::OSIntegral<double,bagel::Int_t::Standard>;
extern template class bagel::OSIntegral<std::complex<double>,bagel::Int_t::London>;
extern template class bagel::OSIntegral<std::complex<double>,bagel::Int_t::Standard>;

#endif

