//
// BAGEL - Parallel electron correlation program.
// Filename: kramers.h
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

#include <src/util/zmatrix.h>
#include <src/util/constants.h>

using namespace std;
using namespace bagel;

namespace bagel {

class Kramers : public ZMatrix {
  protected:

  public:
    Kramers(const int n) : ZMatrix(4*n, 4*n) {

      std::shared_ptr<ZMatrix> unit(new ZMatrix(n, n));
      unit->unit();
      std::shared_ptr<ZMatrix> nunit(new ZMatrix(*unit * -1.0));

      std::complex<double> one  (1.0, 0.0);

      add_block(one, 0, n, n, n, nunit);
      add_block(one, n, 0, n, n, unit);
      add_block(one, 2*n, 3*n, n, n, nunit);
      add_block(one, 3*n, 2*n, n, n, unit);

    }
};

}

