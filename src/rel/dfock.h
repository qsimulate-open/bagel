//
// BAGEL - Parallel electron correlation program.
// Filename: dfock.h
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


#ifndef __SRC_REL_DFOCK_H
#define __SRC_REL_DFOCK_H

#include <memory>
#include <string>
#include <map>
#include <src/wfn/reference.h>
#include <src/wfn/geometry.h>
#include <src/util/zmatrix.h>
#include <src/df/df.h>
#include <src/rel/dfhalfcomplex.h>
#include <src/rel/dfdata.h>
#include <src/rel/relhcore.h>

namespace bagel {

class DFock : public ZMatrix {
  protected:
    std::shared_ptr<const Geometry> geom_;
    const bool gaunt_;

    void two_electron_part(const std::shared_ptr<const ZMatrix> coeff, const bool rhf, const double scale_ex);

    std::list<std::shared_ptr<DFData>> make_dfdists(std::vector<std::shared_ptr<const DFDist>>, bool);
    std::list<std::shared_ptr<DFHalfComplex>> make_half_complex(std::list<std::shared_ptr<DFData>>, std::array<std::shared_ptr<const Matrix>,4>,
                                                                std::array<std::shared_ptr<const Matrix>,4>);

    void add_Jop_block(std::list<std::shared_ptr<DFHalfComplex>>, std::shared_ptr<const DFData>, std::list<std::shared_ptr<const ZMatrix>>);
    void add_Exop_block(std::shared_ptr<DFHalfComplex>, std::shared_ptr<DFHalfComplex>, const double ecale_exch);

  public:
    DFock(const std::shared_ptr<const Geometry> a, 
          const std::shared_ptr<const RelHcore> h,
          const std::shared_ptr<const ZMatrix> coeff, const bool gaunt,
          const bool rhf = true, const double scale_ex = 1.0)
     : ZMatrix(*h), geom_(a), gaunt_(gaunt) {
       
       two_electron_part(coeff, rhf, scale_ex);
    }
    
};

}

#endif
