//
// BAGEL - Parallel electron correlation program.
// Filename: dfhalfbase.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#ifndef __SRC_REL_RELDFBASE_H
#define __SRC_REL_RELDFBASE_H

#include <src/rel/alpha.h>
#include <src/math/zmatrix.h>

namespace bagel {

class SpinorInfo {
  protected:
    // 4*1 matrix
    std::array<std::shared_ptr<const ZMatrix>, 2> spinor_;
    // X and Y (L+, L-, S+, S-) 
    std::pair<int, int> basis_;
    // factor
    std::complex<double> fac_;
    int alpha_comp_;

    void compute_spinor(std::pair<int, int>& cartesian, std::shared_ptr<const Alpha> a) {
      const int index1 = basis_.first;
      const int index2 = basis_.second;
      auto tmp0 = std::make_shared<ZMatrix>(4,1,true);
      auto tmp1 = std::make_shared<ZMatrix>(4,1,true);
      tmp0->element(index1,0) = 1.0;
      tmp1->element(index2,0) = 1.0;
      spinor_[0] = tmp0;
      spinor_[1] = tmp1;

      const Sigma s1(cartesian.first);
      const Sigma s2(cartesian.second);

      ZMatrix z1(s1**spinor_[0]);
      ZMatrix z2(s2**spinor_[1]);
      fac_ = (z1 % *a * z2).element(0,0);
    }
  public:
    SpinorInfo(std::pair<int, int> moblock, std::shared_ptr<const Alpha> alpha, std::pair<int, int> cartesian)
      : basis_(moblock), alpha_comp_(alpha->comp()) {
      compute_spinor(cartesian, alpha);
    }

    SpinorInfo(const SpinorInfo& ab, const int alpha = -1)
      : spinor_(ab.spinor_), basis_(ab.basis_), fac_(ab.fac_), alpha_comp_(alpha == -1 ? ab.alpha_comp_ : alpha) {
    }

    std::complex<double> fac() const { return fac_; }
    bool nonzero() const { return fac_ != std::complex<double>(0.0); }

    bool operator==(const SpinorInfo& o) const { return basis_.first == o.basis_.first && basis_.second == o.basis_.second && fac_ == o.fac_ && alpha_comp_ == o.alpha_comp_; }
    bool operator!=(const SpinorInfo& o) const { return !(*this == o); }

    std::shared_ptr<const SpinorInfo> swap() const {
      auto out = std::make_shared<SpinorInfo>(*this);
      std::swap(out->basis_.first, out->basis_.second);
      out->fac_ = std::conj(out->fac_);
      std::swap(out->spinor_[0], out->spinor_[1]);
      return out;
    }

    int basis(const int i) const {
      assert(i == 0 || i == 1);
      int out = 0;
      for ( ; out != 4; ++out)
        if (spinor_[i]->data(out) == 1.0) break;
      assert(out != 4);
      return out;
    }

    std::pair<int, int> basis() const {return basis_; }
    int basis_first() const { return basis_.first; }
    int basis_second() const { return basis_.second; }
    int comp() const { return alpha_comp_; }
    std::array<std::shared_ptr<const ZMatrix>, 2> spinors() const { return spinor_; }
};


class RelDFBase {
  protected:
    // l,x,y,z
    std::pair<int, int> cartesian_;
    // X,Y, and coefficients
    std::vector<std::shared_ptr<const SpinorInfo>> basis_;

    virtual void set_basis() = 0;

    void common_init() { set_basis(); }

  public:
    RelDFBase(std::pair<int, int> cartesian) : cartesian_(cartesian) {
    }

    RelDFBase(const RelDFBase& o) : cartesian_(o.cartesian_) {
    }

    std::pair<int, int> cartesian() const { return cartesian_; }
    const std::vector<std::shared_ptr<const SpinorInfo>>& basis() const { return basis_; }

};

}

#endif
