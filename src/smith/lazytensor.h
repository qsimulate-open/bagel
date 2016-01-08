//
// BAGEL - Parallel electron correlation program.
// Filename: lazytensor.h
// Copyright (C) 2015 Toru Shiozaki
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


#ifndef __SRC_SMITH_LAZYTENSOR_H
#define __SRC_SMITH_LAZYTENSOR_H

#include <tatensor.h>

namespace bagel {
namespace SMITH {

// derived from the corresponding code in MPQC
// for the time being I assume that LazyTensor is not used with Block expressions
// (i.e., assuming zero-based index range).

template<typename DataType, int N, typename Gen>
class LazyTensor {
  public:
    using value_type = DataType;
    using eval_type = TiledArray::Tensor<DataType>;
    using range_type = typename eval_type::range_type;

    template <typename Archive>
    void serialize(Archive& ar) {
      ar & range_ & gen_;
    }

  private:
    range_type range_;
    Gen gen_;

  public:
    LazyTensor(const range_type& r, const Gen& gen) : range_(r), gen_(gen) {
    }

    LazyTensor(const LazyTensor& other) = default;
    LazyTensor(LazyTensor&& other) = default;
    LazyTensor() = default;

    LazyTensor& operator=(const LazyTensor& other) = default;
    LazyTensor& operator=(LazyTensor&& other) = default;

    operator eval_type () const {
      eval_type tile(range_);
      for (size_t i = 0; i < tile.size(); ++i)
        tile[i] = gen_(i);
      return tile;
    }
};


class GenSample {
  public:
    template <typename Archive>
    void serialize(Archive& ar) {
      // Serialize data members here ...
    }
  protected:
    // members
  public:
    GenSample(} {
    }
    template <typename Index>
    double operator()(const Index& i) {
      return 0.0;
    }
};


template<typename DataType, int N, typename Gen>
class LazyTATensor : public TiledArray::Array<DataType,N,LazyTensor<DataType,N,Gen>> {
  protected:
    using Lazy = LazyTensor<DataType,N,Gen>;
    Gen gen_;

  public:
    using BaseArray = typename TiledArray::Array<DataType,N,Lazy>;
    using BaseArray::begin;
    using BaseArray::end;
    using BaseArray::trange;
    using BaseArray::get_world;

  public:
    LazyTATensor(const std::vector<IndexRange>& r, Gen gen)
      : TiledArray::Array<DataType,N,Lazy>(madness::World::get_default(), *TATensor<DataType,N>::make_trange(r)), gen_(gen) {

      for (auto it = begin(); it != end; ++it) {
        const TiledArray::Range range = trange().make_tile_range(it.ordinal());
        madness::Future<typename Lazy::value_type> tile(Lazy(t, gen_));
        *it = tile;
      }
    }

    auto operator()(const std::string& vars) const -> decltype(TiledArray::Array<DataType,N,Lazy>::operator()(vars)) {
      // due to column-row major covention, we reverse vars
#if defined(__GNUC__) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 9
      using std::regex;
      using std::smatch;
      using std::regex_search;
      using std::sregex_token_iterator;
#else
      using boost::regex;
      using boost::smatch;
      using boost::regex_search;
      using boost::sregex_token_iterator;
#endif
      regex re(",");
      sregex_token_iterator p(vars.begin(), vars.end(), re, -1);
      sregex_token_iterator end;
      std::vector<std::string> ss;
      while (p != end)
        ss.push_back(*p++);
      assert(ss.size() == N);

      std::string revvars;
      auto i = ss.rbegin();
      for (revvars += *i++; i != ss.rend(); ++i)
        revvars += "," + *i;
      return TiledArray::Array<DataType,N,Lazy>::operator()(revvars);
    }
};


}}

#endif
