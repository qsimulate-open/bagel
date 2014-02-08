//
// BAGEL - Parallel electron correlation program.
// Filename: determinants_base.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <iomanip>
#include <src/math/comb.h>
#include <src/fci/fciblock.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::FCIBlockInfo)

using namespace std;
using namespace bagel;

const static Comb comb;

FCIBlockInfo::FCIBlockInfo(const int norb, const int nelea, const int neleb, const bool mute)
  : CIBlockInfo<FCIString>(make_shared<FCIString>(nelea, norb), make_shared<FCIString>(neleb, norb)) {

  if (!mute) {
    cout << "  Performs exactly the same way as Knowles & Handy 1984 CPL" << endl << endl;
    cout << "  o alpha-beta strings" << endl;
    cout << "      length: " << setw(13) << lena() + lenb() << endl;
    cout << "  o size of the space " << endl;
    cout << "      determinant space:  " << lena() * lenb() << endl;
    cout << "      spin-adapted space: " << ncsfs() << endl << endl;
  }
}


size_t FCIBlockInfo::ncsfs() const {
  const int twoS = abs(nspin());
  const int N = nelea() + neleb();
  const int M = norb();
  const size_t out = (twoS + 1) * comb.c(M+1, (N-twoS)/2) * comb.c(M+1, (M-((N+twoS)/2)));
  return out / (M+1);
}
