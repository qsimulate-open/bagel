//
// BAGEL - Parallel electron correlation program.
// Filename: dfock.h
// Copyright (C) 2013 Matthew Kelley
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
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


#ifndef __SRC_REL_DFOCK_H
#define __SRC_REL_DFOCK_H

#include <src/wfn/reference.h>
#include <src/wfn/geometry.h>
#include <src/math/zmatrix.h>
#include <src/rel/reldf.h>
#include <src/rel/relhcore.h>
#include <src/rel/cdmatrix.h>

namespace bagel {

class DFock : public ZMatrix {
  protected:
    std::shared_ptr<const Geometry> geom_;
    const bool gaunt_;
    const bool breit_;

    void two_electron_part(const std::shared_ptr<const ZMatrix> coeff, const bool rhf, const double scale_ex);


    void add_Jop_block(std::shared_ptr<const RelDF>, std::list<std::shared_ptr<const CDMatrix>>, const double scale);
    void add_Exop_block(std::shared_ptr<RelDFHalf>, std::shared_ptr<RelDFHalf>, const double scale, const bool diag = false);
    void driver(std::array<std::shared_ptr<const Matrix>, 4> rocoeff, std::array<std::shared_ptr<const Matrix>, 4> iocoeff,
                           std::array<std::shared_ptr<const Matrix>, 4> trocoeff, std::array<std::shared_ptr<const Matrix>, 4>tiocoeff, bool gaunt, bool breit,
                           const double scale_exchange);

    // when gradient is requested, we store half-transformed integrals
    bool store_half_; 
    std::list<std::shared_ptr<RelDFHalf>> half_; 

  public:
    DFock(const std::shared_ptr<const Geometry> a,
          const std::shared_ptr<const ZMatrix> hc,
          const std::shared_ptr<const ZMatrix> coeff, const bool gaunt, const bool breit,
          const bool store_half, const bool rhf = true, const double scale_ex = 1.0)
     : ZMatrix(*hc), geom_(a), gaunt_(gaunt), breit_(breit), store_half_(store_half) {

       assert(breit ? gaunt : true);
       two_electron_part(coeff, rhf, scale_ex);
    }

    // Utility functions. They are static so that it could be used from gradient codes

    // T needs to have "zaxpy" and "matches" functions.
    template<class T> static void factorize(T& m) {
      for (auto i = m.begin(); i != m.end(); ++i)
        for (auto j = i; j != m.end(); ) {
          if (i != j && (*i)->matches(*j)) {
            (*i)->ax_plus_y((*j)->fac() / (*i)->fac(), *j);
            j = m.erase(j);
          } else
            ++j;
        }
    }
    static std::list<std::shared_ptr<RelDF>> make_dfdists(std::vector<std::shared_ptr<const DFDist>>, bool);
    static std::list<std::shared_ptr<RelDFHalf>> make_half_complex(std::list<std::shared_ptr<RelDF>>, std::array<std::shared_ptr<const Matrix>,4>,
                                                                   std::array<std::shared_ptr<const Matrix>,4>);

    std::list<std::shared_ptr<RelDFHalf>> half() const { assert(store_half_); return half_; }

};

}

#endif
