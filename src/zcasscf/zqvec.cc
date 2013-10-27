//
// BAGEL - Parallel electron correlation program.
// Filename: zqvec.cc
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

#include <src/zcasscf/zqvec.h>

using namespace std;
using namespace bagel;
    

ZQvec::ZQvec(const int n, const int m, shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> coeff, shared_ptr<const ZHarrison> fci)
 : ZMatrix(n,m) {

  array<list<shared_ptr<RelDFHalf>>,2> half_coulomb = fci->jop()->half_complex_coulomb();
  array<list<shared_ptr<RelDFHalf>>,2> half_gaunt   = fci->jop()->half_complex_gaunt();

  array<shared_ptr<const ZMatrix>,2> kcoeff = fci->kramers_coeff();

}
