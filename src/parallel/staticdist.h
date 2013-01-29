//
// BAGEL - Parallel electron correlation program.
// Filename: staticdist.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __SRC_PARALLEL_STATICDIST_H
#define __SRC_PARALLEL_STATICDIST_H

#include <vector>
#include <tuple>
#include <sstream>
#include <stdexcept>

class StaticDist {
  protected:
    const size_t nele_;
    const size_t nproc_;

    std::vector<size_t> start_;

  public:
    StaticDist(const size_t nele, const size_t np) : nele_(nele), nproc_(np) {
      if (nele_ < nproc_) {
        std::stringstream ss;
        ss << "Parallelization with StaticDist is only supported with Nproc smaller than the number of elements. Nele " << nele_ << " Nproc " << nproc_;
        throw std::runtime_error(ss.str());
      }

      const size_t maxsize = (nele_-1) / nproc_ + 1;
      const size_t ares = (nele_-1) % nproc_ + 1;
      for (size_t i = 0; i != nproc_; ++i) {
        start_.push_back((maxsize-1) * i + std::min(i, ares));
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
        out.push_back(std::make_pair(start_[i], start_[i+1]-start_[i]));
      return out;
    }

    std::tuple<size_t, size_t> range(const size_t i) const { assert(i < start_.size()-1); return std::make_tuple(start_[i], start_[i+1]); }
    size_t size(const size_t i) const { return start_[i+1]-start_[i]; }

    std::tuple<size_t, size_t> locate(size_t element) const {
      for (size_t i = 0; i != nproc_; ++i) {
        if (element < start_[i+1]) return std::make_tuple(i, element-start_[i]);
      }
      throw std::runtime_error("wrong call to StaticDist::iproc");
      return std::make_tuple(0,0);
    }

};

#endif
