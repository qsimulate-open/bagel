//
// BAGEL - Parallel electron correlation program.
// Filename: pcoeff.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/periodic/pcoeff.h>

using namespace std;
using namespace bagel;

PCoeff::PCoeff(const int blocksize, const int nblock) : PData(blocksize, nblock) { }

PCoeff::PCoeff(shared_ptr<const Lattice> l) : PData(l->primitive_cell()->nbasis(), l->num_lattice_kvectors()) { }

PCoeff::PCoeff(const PData& o) : PData(o.blocksize(), o.nblock()) {

  pdata_.resize(nblock_);
  for (int i = 0; i != nblock_; ++i)
    pdata_[i] = o.pdata(i);
}

shared_ptr<const PData> PCoeff::form_density_rhf(const int n, const int offset) const {

  PData out(blocksize_, nblock_);
  for (int i = 0; i != nblock_; ++i) {
    const ZMatrix tmp = pdata_[i]->slice(offset, offset + n);
    auto den = make_shared<ZMatrix>(tmp ^ tmp);
    out[i] = den;
  }

  return make_shared<const PData>(out);
}

