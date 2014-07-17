//
// BAGEL - Parallel electron correlation program.
// Filename: product_ci.h
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __BAGEL_ASD_DMRG_PRODUCT_RAS_H
#define __BAGEL_ASD_DMRG_PRODUCT_RAS_H

#include <algorithm>

#include <src/ras/ras_space.h>
#include <src/asd_dmrg/block_key.h>
#include <src/math/matrix.h>
#include <src/ras/civector.h>

namespace bagel {

/// Stores one set of \f$|L\rangle\otimes|\mbox{RAS}\rangle\f$ CI vectors with set combinations of \f$(N_\alpha,N_\beta)_L\otimes(N_\alpha,N_\beta)_{\mbox{RAS}}\f$
class RASBlockVectors : public Matrix {
  protected:
    std::shared_ptr<const RASDeterminants> det_; ///< Determinants information for the \f$|\mbox{RAS}\rangle\f$ part of the vector
    BlockInfo left_state_; ///< /f$|L\rangle\f$

  public:
    /// Constructor
    RASBlockVectors(std::shared_ptr<const RASDeterminants> det, const BlockInfo ls) : Matrix(det->size(), ls.nstates), det_(det), left_state_(ls) {
      assert(det);
    }
    /// Copy-constructor
    RASBlockVectors(const RASBlockVectors& o) : RASBlockVectors(o.det_, o.left_state_) {
      std::copy_n(o.data(), size(), data());
    }

    std::vector<RASCivecView> civecs() {
      std::vector<RASCivecView> out;
      const int nst = nstates();
      for (int i = 0; i < nst; ++i)
        out.emplace_back(det_, element_ptr(0,i));
      return out;
    }
    RASCivecView civec(const int i) { return RASCivecView(det_, element_ptr(0,i)); }

    int nstates() const { return left_state_.nstates; }

    using Matrix::data;
    using Matrix::element_ptr;
};

/**
  ProductRASCivec stores a CI vector of the type \f$|L\rangle\otimes|\mbox{RAS}\rangle\f$ where \f$|L\rangle\f$ is specified with BlockInfo
  and \f$|\mbox{RAS}\rangle\f$ is a regular RAS wavefunction
*/
class ProductRASCivec {
  protected:
    std::map<BlockKey, std::shared_ptr<RASBlockVectors>> sectors_; ///< The key for sectors_ specifies information for the left_block
    std::shared_ptr<RASSpace> space_;

    std::set<BlockInfo> lblocks_;

    int nelea_;
    int neleb_;

  public:
    ProductRASCivec(std::shared_ptr<RASSpace> space, std::set<BlockInfo>& left_blocks, const int nelea, const int neleb); ///< Constructor
    ProductRASCivec(const ProductRASCivec& o); ///< Copy-constructor
    ProductRASCivec(ProductRASCivec&& o); ///< Move-constructor

    ProductRASCivec& operator=(const ProductRASCivec& o); ///< Copy-assignment
    ProductRASCivec& operator=(ProductRASCivec&& o); ///< Move-assignment

    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }

    size_t size() const {
      return std::accumulate(sectors_.begin(), sectors_.end(), 0ul, [] (const size_t a, const std::pair<const BlockKey, std::shared_ptr<RASBlockVectors>>& p) { return a+p.second->size(); });
    }

    const std::set<BlockInfo>& lblocks() const { return lblocks_; }

    bool matches(const ProductRASCivec& o) const {
      return (*space_==*o.space_ && std::make_pair(nelea_,neleb_)==std::make_pair(o.nelea(),o.neleb()) && std::equal(lblocks_.begin(), lblocks_.end(), o.lblocks_.begin()));
    }

    void scale(const double a);
    void ax_plus_y(const double& a, const ProductRASCivec& o);

    double dot_product(const ProductRASCivec& o) const ;

    double norm() const { return std::sqrt(dot_product(*this)); }
    double variance() const { return dot_product(*this)/size(); }
    double rms() const { return std::sqrt(variance()); }
};

}

#endif
