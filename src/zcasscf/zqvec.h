//
// BAGEL - Parallel electron correlation program.
// Filename: zqvec.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __BAGEL_SRC_ZCASSCF_ZQVEC_H
#define __BAGEL_SRC_ZCASSCF_ZQVEC_H

#include <src/zfci/zharrison.h> // 2RDM and transformed integrals

namespace bagel {

class ZQvec : public ZMatrix {
  protected:

  public:
    ZQvec(const int n, const int m, std::shared_ptr<const Geometry> geom, std::shared_ptr<const ZMatrix> coeff, const int nclosed,
          std::shared_ptr<const ZHarrison> fci, const bool gaunt, const bool breit); // FIXME : this constructor has a bug

    ZQvec(const ZMatrix& a) : ZMatrix(a) {}

    ZQvec(const int n, const int m, std::shared_ptr<const Geometry> geom, std::shared_ptr<const ZMatrix> rcoeff, std::shared_ptr<const ZMatrix> acoeff, const int nclosed,
          std::shared_ptr<const ZHarrison> fci, const bool gaunt, const bool breit);
};

}

#endif
