//
// BAGEL - Parallel electron correlation program.
// Filename: loopgenerator.cc
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/loopgenerator.h>
#include <algorithm>

using namespace bagel::SMITH;
using namespace std;


vector<vector<Index>> LoopGenerator::block_loop() const {
  // first, make a status vector
  vector<int> stat(loop_.size());
  vector<int> max(loop_.size());
  {
    auto j = max.begin();
    for (auto& i : loop_) *j++ = i.nblock();
  }

  vector<vector<Index>> out;

  do {
    vector<Index> tmp;
    auto k = stat.begin();
    for (auto& l : loop_)
      tmp.push_back(l.range(*k++));
    out.push_back(tmp);

    auto j = stat.begin();
    auto i = max.begin();
    while (j != stat.end() && (++*j) == *i) {
      *j = 0;
      ++j;
      ++i;
    }
  } while (*max_element(stat.begin(), stat.end()) > 0);

  return out;
}

#endif
