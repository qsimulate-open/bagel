//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: product_citask.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __BAGEL_ASD_DMRG_PRODUCT_CITASK_H
#define __BAGEL_ASD_DMRG_PRODUCT_CITASK_H

namespace bagel {

/// Base for computing columns of model operators in the ProductCI framework
template <typename Derived>
class ProductCITask {
  protected:
    std::vector<PCI::Basis>* basis_;
    const size_t col1_;
    const size_t col2_;
    double* const dest1_;
    double* const dest2_;


    ProductCITask(std::vector<PCI::Basis>* b, const size_t c1, double* d1, const size_t c2, double* d2) :
      basis_(b), col1_(c1), col2_(c2), dest1_(d1), dest2_(d2) { }

    double matrix_element(const PCI::Basis& bra, const PCI::Basis& ket) { return static_cast<Derived*>(this)->matrix_element_impl(bra, ket); }

  public:
    void compute() {
      double* odata = dest1_;
      const PCI::Basis bra1 = basis_->at(col1_);
      std::for_each(basis_->begin() + col1_, basis_->end(), [this, &bra1, &odata] (const PCI::Basis& ket) { *odata = matrix_element(bra1, ket); ++odata; });

      if (col1_ != col2_) {
        odata = dest2_;
        const PCI::Basis bra2 = basis_->at(col2_);
        std::for_each(basis_->begin() + col2_, basis_->end(), [this, &bra2, &odata] (const PCI::Basis& ket) { *odata = matrix_element(bra2, ket); ++odata; });
      }
    }
};

class ProductCIHamTask : public ProductCITask<ProductCIHamTask> {
  protected:
    std::shared_ptr<const BlockOperators> blockops_;
    std::shared_ptr<const DimerJop> jop_;
    std::shared_ptr<const Matrix> mo1e_;

    const int rnorb_;

    double mo2e(int i, int j, int k, int l) {
      if (i > j) std::swap(i,j);
      if (k > l) std::swap(k,l);
      return jop_->mo2e(i,j,k,l);
    }

    double compute_pure_ras(const std::bitset<nbit__> abra, const std::bitset<nbit__> bbra, const std::bitset<nbit__> aket, const std::bitset<nbit__> bket);

  public:
    ProductCIHamTask(std::vector<PCI::Basis>* b, std::shared_ptr<const BlockOperators> blockops, std::shared_ptr<const DimerJop> jop, std::shared_ptr<const Matrix> mo1e, const size_t c1, double* d1, const size_t c2, double* d2);

    double matrix_element_impl(PCI::Basis bra, PCI::Basis ket);
};

}

#endif
