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
#include <src/util/zmatrix.h>

namespace bagel {

class ABcases {
  protected:
    std::array<std::shared_ptr<const ZMatrix>, 2> spinor_;
    std::pair<int, int> basis_; 
    std::complex<double> fac_;
    int alpha_comp_;

    void compute_spinor(std::pair<int, int>& coord, std::shared_ptr<const Sigma> s1, std::shared_ptr<const Sigma> s2, std::shared_ptr<const Alpha> a) {
      const int start1 = coord.first == Comp::L ? 0 : 2;
      const int start2 = coord.second == Comp::L ? 0 : 2;
      const int index1 = start1 + basis_.first;
      const int index2 = start2 + basis_.second;
      std::shared_ptr<ZMatrix> tmp0(new ZMatrix(4,1,true));
      std::shared_ptr<ZMatrix> tmp1(new ZMatrix(4,1,true));
      tmp0->element(index1,0) = 1.0;
      tmp1->element(index2,0) = 1.0;
      spinor_[0] = tmp0;
      spinor_[1] = tmp1;

      ZMatrix z1(*s1->data()**spinor_[0]);
      ZMatrix z2(*s2->data()**spinor_[1]);
      fac_ = (z1 % *a->data() * z2).element(0,0);
    }
  public:
    ABcases(std::pair<int, int> b, std::pair<int, int> c, std::shared_ptr<const Sigma> s1, std::shared_ptr<const Sigma> s2, std::shared_ptr<const Alpha> a)
      : basis_(b), alpha_comp_(a->comp()) {
      compute_spinor(c, s1, s2, a);
    }

    ABcases(const ABcases& ab, const int alpha = -1)
      : spinor_(ab.spinor_), basis_(ab.basis_), fac_(ab.fac_), alpha_comp_(alpha == -1 ? ab.alpha_comp_ : alpha) {
    }

    std::complex<double> fac() const { return fac_; }
    bool nonzero() const { return fac_ != std::complex<double>(0.0); }

    bool operator==(const ABcases& o) const { return basis_.first == o.basis_.first && basis_.second == o.basis_.second && fac_ == o.fac_ && alpha_comp_ == o.alpha_comp_; } 
    bool operator!=(const ABcases& o) const { return !(*this == o); }

    std::shared_ptr<const ABcases> swap() const {
      std::shared_ptr<ABcases> out(new ABcases(*this));
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
    std::pair<int, int> coord_;
    std::vector<std::shared_ptr<const ABcases>> basis_;

    virtual void set_basis() = 0;

    void common_init() { set_basis(); }

  public:
    RelDFBase(std::pair<int, int> coord) : coord_(coord) {
    }

    RelDFBase(const RelDFBase& o) : coord_(o.coord_) {
    }

    std::pair<int, int> coord() const { return coord_; }
    const std::vector<std::shared_ptr<const ABcases>>& basis() const { return basis_; }

};

}

#endif
