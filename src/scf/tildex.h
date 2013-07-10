//
// BAGEL - Parallel electron correlation program.
// Filename: tildex.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#ifndef __src_scf_tildex_h
#define __src_scf_tildex_h

#include <src/math/matrix.h>

namespace bagel {

class TildeX : public Matrix {
  protected:

  public:
    TildeX(const std::shared_ptr<const Overlap> olp, const double thresh) : Matrix(*olp)  {
      this->inverse_half(thresh);
    }

    TildeX& operator=(const Matrix1e& o) { std::copy(o.data(), o.data()+ndim_*mdim_, this->data()); return *this; };

};

}

#endif
