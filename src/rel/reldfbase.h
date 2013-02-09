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

class RelDFBase {
  protected:
    std::pair<int, int> coord_;
    std::pair<int, int> basis_; 
    std::shared_ptr<const Alpha> alpha_;
    std::shared_ptr<const Sigma> sigma1_;
    std::shared_ptr<const Sigma> sigma2_;
    std::array<std::shared_ptr<ZMatrix>, 2> spinor_;
    std::complex<double> fac_;

    void compute_spinor() {
      const int start1 = coord_.first == Comp::L ? 0 : 2;
      const int start2 = coord_.second == Comp::L ? 0 : 2;
      const int index1 = start1 + basis_.first;
      const int index2 = start2 + basis_.second;
      spinor_[0] = std::shared_ptr<ZMatrix>(new ZMatrix(4,1,true));
      spinor_[1] = std::shared_ptr<ZMatrix>(new ZMatrix(4,1,true));
      spinor_[0]->element(index1,0) = 1.0;
      spinor_[1]->element(index2,0) = 1.0;
    }

    void common_init() {
      sigma1_ = std::shared_ptr<const Sigma>(new Sigma(coord_.first));
      sigma2_ = std::shared_ptr<const Sigma>(new Sigma(coord_.second));

      compute_spinor();

      ZMatrix z1(*sigma1_->data()**spinor_[0]);
      ZMatrix z2(*sigma2_->data()**spinor_[1]);
      fac_ = (z1 % *alpha_->data() * z2).element(0,0);

      assert(fac_ != std::complex<double>(0.0));
    }

  public:
    RelDFBase(std::pair<int, int> coord, const int alpha) : coord_(coord) {
      if ((coord_.first == Comp::Z) ^ (coord_.second == Comp::Z))
        basis_ = std::make_pair(Basis::a, Basis::b);
      else
        basis_ = std::make_pair(Basis::a, Basis::a);

      alpha_ = std::shared_ptr<const Alpha>(new Alpha(alpha));
      common_init();
    }

    RelDFBase(const RelDFBase& o) : coord_(o.coord_), basis_(o.basis_), alpha_(o.alpha_) {
      common_init();
    }


    std::pair<int, int> coord() const { return coord_; }
    std::pair<int, int> basis() const { return basis_; }

    std::shared_ptr<const Alpha> alpha() const { return alpha_; }
    std::complex<double> fac() const { return fac_; }

    int basis(const int i) const {
      assert(i == 0 || i == 1);
      int out = 0;
      for ( ; out != 4; ++out)
        if (spinor_[i]->data(out) == 1.0) break;
      assert(out != 4);
      return out;
    }

};

}

#endif
