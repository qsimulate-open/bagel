//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: indexrange.cc
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

#include <src/smith/indexrange.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


IndexRange::IndexRange(const int size, const int maxblock, const int boffset, const int orboff, const int orboff2)
 : keyoffset_(boffset), orboffset_(orboff), orboffset2_(orboff2 < 0 ? orboff : orboff2) {
  if (size > 0) {
    // first determine number of blocks.
    const size_t nbl = (size-1) / maxblock + 1;
    const size_t nblock = (size-1) / nbl + 1;
    // we want to distribute orbitals as evenly as possible
    const size_t rem = nbl * nblock - size;
    vector<size_t> blocksizes(nbl, nblock);
    auto iter = blocksizes.rbegin();
    for (int k = 0; k != rem; ++iter, ++k) --*iter;
    // push back to range_
    size_t off = orboffset_;
    size_t off2 = orboffset2_;
    // key is offsetted by the input value
    size_t cnt = boffset;
    for (auto& i : blocksizes) {
      range_.emplace_back(off, off2, i, cnt++);
      off += i;
      off2 += i;
    }
    // set size_
    size_ = off - orboffset_;
  } else {
    size_ = 0;
  }
}


IndexRange::IndexRange(shared_ptr<const StaticDist> dist, const int boffset, const int orboff, const int orboff2)
 : keyoffset_(boffset), orboffset_(orboff), orboffset2_(orboff2 < 0 ? orboff : orboff2) {
  // push back to range_
  size_t off = orboffset_;
  size_t off2 = orboffset2_;
  // key is offsetted by the input value
  size_t cnt = boffset;
  vector<pair<size_t,size_t>> table = dist->atable();
  for (auto& i : table) {
    range_.emplace_back(off, off2, i.second, cnt++);
    off += i.second;
    off2 += i.second;
  }
  // set size_
  size_ = off - orboffset_;
}


void IndexRange::merge(const IndexRange& o) {
  range_.insert(range_.end(), o.range_.begin(), o.range_.end());
  size_ += o.size_;
}


bool IndexRange::operator==(const IndexRange& o) const {
  bool out = size_ == o.size_;
  if (range_.size() == o.range_.size()) {
    auto i = range_.begin();
    for (auto j = o.range_.begin(); i != range_.end(); ++i, ++j) out &= (*i) == (*j);
  } else {
    out = false;
  }
  return out;
}


string IndexRange::str() const {
  stringstream ss;
  for (auto& i : range_)
    ss << setw(10) << i.offset() << setw(10) << i.size() << endl;
  return ss.str();
}

BOOST_CLASS_EXPORT_IMPLEMENT(Index)
BOOST_CLASS_EXPORT_IMPLEMENT(IndexRange)
