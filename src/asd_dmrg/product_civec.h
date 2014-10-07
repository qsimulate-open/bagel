//
// BAGEL - Parallel electron correlation program.
// Filename: product_civec.h
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
#include <src/asd_dmrg/dmrg_block.h>
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

    /// Returns a set of RASCivecView objects that interpret the BlockVector as CI vectors
    std::vector<RASCivecView> civecs() {
      std::vector<RASCivecView> out;
      const int nst = nstates();
      for (int i = 0; i < nst; ++i)
        out.emplace_back(det_, element_ptr(0,i));
      return out;
    }
    /// Interpret the i-th vector as a RASCivecView
    RASCivecView civec(const int i) { return RASCivecView(det_, element_ptr(0,i)); }
    /// Interpret the i-th vector as a const RASCivecView
    const RASCivecView civec(const int i) const { return RASCivecView(det_, element_ptr(0,i)); }

    BlockInfo left_state() const { return left_state_; }
    int nstates() const { return left_state_.nstates; }
    std::shared_ptr<const RASDeterminants> det() const { return det_; }

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
    std::shared_ptr<const DMRG_Block> left_;

    int nelea_;
    int neleb_;

  public:
    ProductRASCivec(std::shared_ptr<RASSpace> space, std::shared_ptr<const DMRG_Block> left, const int nelea, const int neleb); ///< Constructor
    ProductRASCivec(const ProductRASCivec& o); ///< Copy-constructor
    ProductRASCivec(ProductRASCivec&& o); ///< Move-constructor

    ProductRASCivec& operator=(const ProductRASCivec& o); ///< Copy-assignment
    ProductRASCivec& operator=(ProductRASCivec&& o); ///< Move-assignment

    std::shared_ptr<ProductRASCivec> clone() const { return std::make_shared<ProductRASCivec>(space_, left_, nelea_, neleb_); }
    std::shared_ptr<ProductRASCivec> copy() const { return std::make_shared<ProductRASCivec>(*this); }

    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
    int nele()  const { return nelea_ + neleb_; }

    /// Total size of configuration space.
    size_t size() const {
      return std::accumulate(sectors_.begin(), sectors_.end(), 0ul, [] (const size_t a, const std::pair<const BlockKey, std::shared_ptr<RASBlockVectors>>& p) { return a+p.second->size(); });
    }

    std::shared_ptr<const DMRG_Block> left() const { return left_; }
    std::shared_ptr<RASSpace> space() { return space_; }
    std::shared_ptr<const RASSpace> space() const { return space_; }

    /// Check whether a block with the given key exists.
    bool contains_block(const BlockKey bk) const { return sectors_.find(bk)!=sectors_.end(); }

    std::map<BlockKey, std::shared_ptr<RASBlockVectors>>& sectors() { return sectors_; }
    const std::map<BlockKey, std::shared_ptr<RASBlockVectors>>& sectors() const { return sectors_; }

    std::shared_ptr<RASBlockVectors> sector(const BlockKey& b) { return sectors_.at(b); }
    std::shared_ptr<const RASBlockVectors> sector(const BlockKey& b) const { return sectors_.at(b); }

    /// Check whether two ProductRASCivecs are compatible.
    /** ProductRASCivecs are compatible if they have the same RAS structure,
        the same total number of electrons, and the same DMRG_Block. */
    bool matches(const ProductRASCivec& o) const {
      return (*space_==*o.space_ && std::make_pair(nelea_,neleb_)==std::make_pair(o.nelea(),o.neleb()) && left_->blocks()==o.left_->blocks());
    }

    double spin_expectation() const;
    void spin_decontaminate(const double thresh = 1.0e-8);

    std::shared_ptr<ProductRASCivec> spin() const;
    std::shared_ptr<ProductRASCivec> spin_lower() const;
    std::shared_ptr<ProductRASCivec> spin_raise() const;

    void scale(const double a);
    void ax_plus_y(const double& a, const ProductRASCivec& o);
    void ax_plus_y(const double& a, std::shared_ptr<const ProductRASCivec>& o) { ax_plus_y(a, *o); }

    double dot_product(const ProductRASCivec& o) const;
    double dot_product(std::shared_ptr<const ProductRASCivec>& o) const { return dot_product(*o); }

    double norm() const { return std::sqrt(dot_product(*this)); }
    double variance() const { return dot_product(*this)/size(); }
    double rms() const { return std::sqrt(variance()); }

    double normalize() {
      const double nrm = norm();
      scale(nrm > numerical_zero__ ? 1.0/nrm : 0.0);
      return nrm;
    }

    void allreduce();

    void print(const double thresh = 0.05) const;

    /// Splits a vector by the "right" block.
    /** Only defined if left_ is a DMRG_Block2 object, in which case it splits the vector
        into a map such that the right block key is the map key and the mapped object
        is a \f$|L\rangle \otimes |\mbox{RAS}\rangle\f$ type vector. */
    std::map<BlockKey, std::vector<std::shared_ptr<ProductRASCivec>>> split() const;
};

}

#endif
