//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: staticdist.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __SRC_PARALLEL_STATICDIST_H
#define __SRC_PARALLEL_STATICDIST_H

#include <vector>
#include <tuple>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cassert>

namespace bagel {

class StaticDist {
  protected:
    const size_t nele_;
    const size_t nproc_;

    std::vector<size_t> start_;

  public:
    StaticDist(const size_t nele, const size_t np, const size_t chunk = 1) : nele_(nele), nproc_(np) {
      assert(nele % chunk == 0 && nele / chunk > 0);
      const size_t ne = nele_ / chunk;
      const size_t maxsize = (ne-1) / nproc_ + 1;
      const size_t ares = (ne-1) % nproc_ + 1;
      for (size_t i = 0; i != nproc_; ++i) {
        start_.push_back(((maxsize-1) * i + std::min(i, ares)) * chunk);
      }
      start_.push_back(nele);
    }

    // a vector of start and size
    StaticDist(const std::vector<size_t>& o) : nele_(o.back()), nproc_(o.size()-1), start_(o) { assert(o.size() > 1); }

    StaticDist() = delete;

    // vector of pairs of astart and asize
    std::vector<std::pair<size_t, size_t>> atable() const {
      std::vector<std::pair<size_t, size_t>> out;
      for (size_t i = 0; i != nproc_; ++i)
        out.push_back({start_[i], start_[i+1]-start_[i]});
      return out;
    }

    std::tuple<size_t, size_t> range(const size_t i) const { assert(i < start_.size()-1); return std::make_tuple(start_[i], start_[i+1]); }
    size_t start(const size_t i) const { return start_[i]; }
    size_t size(const size_t i) const { return start_[i+1]-start_[i]; }

    size_t nele() const { return nele_; }

    std::tuple<size_t, size_t> locate(size_t element) const {
      for (size_t i = 0; i != nproc_; ++i) {
        if (element < start_[i+1]) return std::make_tuple(i, element-start_[i]);
      }
      throw std::runtime_error("wrong call to StaticDist::iproc");
      return std::make_tuple(0,0);
    }

    void print() const {
      std::cout << "  === StaticDist Information ===" << std::endl;
      for (int i = 0; i != nproc_; ++i)
        std::cout << "   " << std::setw(5) << i << std::setw(8) << start_[i] << "--"
                                                << std::setw(8) << start_[i+1]-1 << std::endl;
    }

};

}

#endif
