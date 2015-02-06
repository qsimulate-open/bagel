//
// BAGEL - Parallel electron correlation program.
// Filename: multi/rasscf/qvec.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <shiozaki@northwestern.edu>
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


#ifndef __BAGEL_SRC_RASSCF_QVEC
#define __BAGEL_SRC_RASSCF_QVEC

#include <src/ci/ras/rasci.h> // 2RDM and half-transformed integrals
#include <src/multi/rasscf/rotfile.h>

namespace bagel {

class RASQvec : public Matrix {
  protected:

  public:
    RASQvec(const int n, const int m, std::shared_ptr<const Matrix> c, const size_t nclosed,
         std::shared_ptr<const RASCI> fci, std::shared_ptr<const RDM<2>> rdm);
    RASQvec(const Matrix& a) : Matrix(a) {}

};

}

#endif
