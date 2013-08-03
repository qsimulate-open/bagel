//
// BAGEL - Parallel electron correlation program.
// Filename: alpha.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Matthew Kelley <matthewkelley2017@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#ifndef __SRC_REL_ALPHA_H
#define __SRC_REL_ALPHA_H

#include <src/math/zmatrix.h>
#include <src/util/constants.h>

namespace bagel {

// to make readable
enum Comp { X = 0, Y = 1, Z = 2, L = 3 };
enum Basis { LP = 0, LM = 1, SP = 2, SM = 3}; 

class Sigma2 : public ZMatrix {
  protected:
  public:
    Sigma2(const int i) : ZMatrix(2,2,true) {
      if (i == Comp::L) {
        element(0,0) = element(1,1) = 1.0;
      } else if (i == Comp::X) {
        element(0,1) = element(1,0) = 1.0;
      } else if (i == Comp::Y) {
        element(1,0) = std::complex<double>(0.0, 1.0);
        element(0,1) = std::complex<double>(0.0,-1.0);
      } else if (i == Comp::Z) {
        element(0,0) = 1.0;
        element(1,1) = -1.0;
      } else {
        assert(false);
      }
    }
};


class Sigma : public ZMatrix {
  protected:
  public:
    Sigma(const int i) : ZMatrix(4,4,true) {
      auto s = std::make_shared<Sigma2>(i);
      if (i == Comp::L) {
        copy_block(0,0,2,2,s);
      } else {
        copy_block(2,2,2,2,s);
        scale(std::complex<double>(0.0, -0.5/c__));
      }
    }
};

class Alpha : public ZMatrix {
  protected:
    const int alpha_comp_;
  public:
    Alpha(const int i) : ZMatrix(4,4,true), alpha_comp_(i) {
      auto s = std::make_shared<Sigma2>(i);
      if (i == Comp::L) {
        copy_block(0, 0, 2, 2, s);
        copy_block(2, 2, 2, 2, s);
      } else if (i == Comp::Z) {
        copy_block(2, 0, 2, 2, s);
        copy_block(0, 2, 2, 2, s);
      } else if (i == Comp::X) {
        copy_block(2, 0, 2, 2, s);
        copy_block(0, 2, 2, 2, s);
      } else if (i == Comp::Y) {
        copy_block(2, 0, 2, 2, s);
        copy_block(0, 2, 2, 2, s);
      } else {
        assert(false);
      }
    }

    const int comp() const { return alpha_comp_; }

};

}

#endif
