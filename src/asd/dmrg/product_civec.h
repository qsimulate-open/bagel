//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: product_civec.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifndef __BAGEL_ASD_DMRG_PRODUCT_RAS_H
#define __BAGEL_ASD_DMRG_PRODUCT_RAS_H

#include <algorithm>
#include <src/ci/ras/ras_space.h>
#include <src/ci/ras/civector.h>
#include <src/asd/dmrg/dmrg_block.h>
#include <src/util/math/matrix.h>

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
    RASBlockVectors(std::shared_ptr<const RASDeterminants> det, const int M) : Matrix(det->size(), M), det_(det), left_state_(0, 0, M) {
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

    RASBlockVectors transpose_civecs(std::shared_ptr<const RASDeterminants> transdet = nullptr) const;

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
      return std::accumulate(sectors_.begin(), sectors_.end(), 0ull, [] (const size_t a, const std::pair<const BlockKey, std::shared_ptr<RASBlockVectors>>& p) { return a+p.second->size(); });
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
    void broadcast(int rank = 0);
    void synchronize();

    void print(const double thresh = 0.05) const;

    /// Splits a vector by the "right" block.
    /** Only defined if left_ is a DMRG_Block2 object, in which case it splits the vector
        into a map such that the right block key is the map key and the mapped object
        is a \f$|L\rangle \otimes |\mbox{RAS}\rangle\f$ type vector. */
    std::map<BlockKey, std::vector<std::shared_ptr<ProductRASCivec>>> split() const;
};

}

#endif
