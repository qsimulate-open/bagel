//
// BAGEL - Parallel electron correlation program.
// Filename: pdata.cc
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

#include <src/periodic/pdata.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(PData)

PData::PData(const int bsize, const int nblock) : blocksize_(bsize), nblock_(nblock) {

  pdata_.resize(nblock);
  for (int i = 0; i != nblock; ++i) {
    auto block = make_shared<ZMatrix>(bsize, bsize);
    block->zero();
    pdata_[i] = block;
  }
}

PData::PData(const PData& o) : blocksize_(o.blocksize()), nblock_(o.nblock()) {

  pdata_.resize(nblock_);
  for (int i = 0; i != nblock_; ++i)
    pdata_[i] = o.pdata(i);

}

shared_ptr<const PData> PData::ft(const vector<array<double, 3>> gvector, const vector<array<double, 3>> kvector) const {

  PData out(blocksize_, kvector.size());

  int k = 0;
  for (auto& kvec : kvector) {
    auto kblock = make_shared<ZMatrix>(blocksize_, blocksize_);
    kblock->zero();
    int g = 0;
    for (auto& gvec : gvector) {
      shared_ptr<const ZMatrix> gblock = pdata(g);
      complex<double> factor(0.0, gvec[0]* kvec[0] + gvec[1] * kvec[1] + gvec[2] * kvec[2]);
      factor = exp(factor);
      auto tmp1 = make_shared<ZMatrix>(*(gblock->get_real_part()), factor); // gblock should be real
      auto tmp2 = make_shared<ZMatrix>(*(gblock->get_imag_part()), factor);
      *kblock += *tmp1 + *tmp2;
      ++g;
    }
    out[k] = kblock;
    ++k;
  }

  return make_shared<const PData>(out);
}

shared_ptr<const PData> PData::ift(const vector<array<double, 3>> gvector, const vector<array<double, 3>> kvector) const {

  PData out(blocksize_, gvector.size());

  int g = 0;
  for (auto& gvec : gvector) {
    auto gblock = make_shared<ZMatrix>(blocksize_, blocksize_);
    gblock->zero();
    int k = 0;
    for (auto& kvec : kvector) {
      shared_ptr<const ZMatrix> kblock = pdata(k);
      complex<double> factor(0.0, -gvec[0]* kvec[0] - gvec[1] * kvec[1] - gvec[2] * kvec[2]);
      factor = exp(factor);
      auto tmp1 = make_shared<ZMatrix>(*(kblock->get_real_part()), factor);
      auto tmp2 = make_shared<ZMatrix>(*(kblock->get_imag_part()), factor);
      *gblock += *tmp1 + *tmp2;
      ++k;
    }
    out[g] = gblock;
    ++g;
  }

  return make_shared<const PData>(out);
}
