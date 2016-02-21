//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: loopgenerator.h
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


#ifndef __SRC_SMITH_LOOPGENERATOR_H
#define __SRC_SMITH_LOOPGENERATOR_H

#include <vector>
#include <src/smith/indexrange.h>

namespace bagel {
namespace SMITH {

struct LoopGenerator {
  static std::vector<std::vector<Index>> gen(const std::vector<IndexRange>& loop) {
    std::vector<int> max(loop.size());
    auto m = max.begin();
    for (auto& i : loop)
      *m++ = i.nblock();

    std::vector<std::vector<Index>> out;
    std::vector<int> stat(loop.size());

    do {
      std::vector<Index> tmp;
      auto k = stat.begin();
      for (auto& l : loop)
        tmp.push_back(l.range(*k++));
      out.push_back(tmp);

      auto j = stat.begin();
      auto i = max.begin();
      while (j != stat.end() && (++*j) == *i) {
        *j = 0;
        ++j;
        ++i;
      }
    } while (*std::max_element(stat.begin(), stat.end()) > 0);

    return out;
  }
};

}
}

#endif
