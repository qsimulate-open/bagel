//
// BAGEL - Parallel electron correlation program.
// Filename: qvec.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __BAGEL_SRC_CASSCF_QVEC
#define __BAGEL_SRC_CASSCF_QVEC

#include <src/fci/fci.h> // 2RDM and half-transformed integrals
#include <src/casscf/rotfile.h>

namespace bagel {

class Qvec : public Matrix {
  protected:

  public:
    Qvec(const int n, const int m, std::shared_ptr<const Matrix> c, const size_t nclosed,
         std::shared_ptr<const FCI> fci, std::shared_ptr<const RDM<2>> rdm);
    Qvec(const Matrix& a) : Matrix(a) {}

};

}

#endif
