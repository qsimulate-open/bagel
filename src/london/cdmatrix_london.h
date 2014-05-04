//
// BAGEL - Parallel electron correlation program.
// Filename: cdmatrix_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_LONDON_CDMATRIX_LONDON_H
#define __SRC_LONDON_CDMATRIX_LONDON_H

#include <src/london/reldfhalf_london.h>
#include <src/math/zmatrix.h>

namespace bagel {

class RelDFHalf_London;

class CDMatrix_London : public ZMatrix {
  protected:
    const int alpha_comp_;

  public:
    CDMatrix_London(std::shared_ptr<const RelDFHalf_London> dfhc, std::shared_ptr<const SpinorInfo> abc, std::array<std::shared_ptr<const ZMatrix>, 4> trcoeff,
                    std::array<std::shared_ptr<const ZMatrix>, 4> ticoeff, std::shared_ptr<const ZMatrix> dat2, const bool onlyonce = true);
    CDMatrix_London(const ZMatrix& o, const int acomp) : ZMatrix(o), alpha_comp_(acomp) { }

    int alpha_comp() const { return alpha_comp_; }

};

}

#endif
