//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: qvec.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __BAGEL_SRC_CASSCF_QVEC
#define __BAGEL_SRC_CASSCF_QVEC

#include <src/ci/fci/fci.h> // 2RDM and half-transformed integrals
#include <src/multi/casscf/rotfile.h>

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
