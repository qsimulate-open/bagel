//
// BAGEL - Parallel electron correlation program.
// Filename: cdmatrix.h
// Copyright (C) 2013 Matthew Kelley
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
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


#ifndef __SRC_REL_CDMATRIX_H
#define __SRC_REL_CDMATRIX_H

#include <src/util/zmatrix.h>

namespace bagel {

class CDMatrix : public ZMatrix {
  protected:
    const int comp_;

  public:
    CDMatrix(const ZMatrix& o, const int comp) : ZMatrix(o), comp_(comp) { }
    const int comp() const { return comp_; }

};

}

#endif
