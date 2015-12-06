//
// BAGEL - Parallel electron correlation program.
// Filename: tatensor.h
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

#ifndef __SRC_SMITH_TATENSOR_H
#define __SRC_SMITH_TATENSOR_H

#include <tiledarray.h>
#include <src/smith/indexrange.h>

#if defined(__GNUC__) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 9
#include <regex>
#else
#include <boost/regex.hpp>
#endif

namespace TiledArray {

// modified version of clone function (which is copy in BAGEL's language)
// see the original in tiledarray/src/TiledArray/conversions/clone.h written by Justus Calvin
template <typename Tile, typename Policy>
inline DistArray<Tile, Policy> clone_part(const DistArray<Tile, Policy>& arg) {
  using value_type = typename DistArray<Tile, Policy>::value_type;
  World& world = arg.get_world();
  // TODO until sparsity is handled properly
  world.gop.fence();
  DistArray<Tile, Policy> result(world, arg.trange(), arg.get_shape(), arg.get_pmap());

  for (auto index : *arg.get_pmap()) {
    if (arg.is_zero(index))
      continue;
    auto atile = arg.find(index);
    if (!atile.probe())
      continue;
    Future<value_type> tile = world.taskq.add([](const value_type& tile){ return TiledArray::clone(tile); }, atile);
    result.set(index, tile);
  }
  // TODO until sparsity is handled properly
  world.gop.fence();
  return result;
}

}

namespace bagel {
namespace SMITH {

// wrapper class for TiledArray Tensor
template<typename DataType, int N>
class TATensor : public TiledArray::Array<DataType,N> {
  public:
    using BaseArray = typename TiledArray::Array<DataType,N>;
    using value_type = typename BaseArray::value_type;
    using size_type  = typename BaseArray::size_type;

    using BaseArray::begin;
    using BaseArray::end;
    using BaseArray::trange;
    using BaseArray::get_world;
    using BaseArray::find;
    using BaseArray::id;

    static std::shared_ptr<TiledArray::TiledRange> make_trange(const std::vector<IndexRange>& r) {
      std::vector<TiledArray::TiledRange1> ranges;
      for (auto it = r.rbegin(); it != r.rend(); ++it) {
        std::vector<size_t> tile_boundaries;
        size_t off = 0;
        for (auto& j : *it) {
          tile_boundaries.push_back(off);
          off += j.size();
        }
        tile_boundaries.push_back(off);
        ranges.emplace_back(tile_boundaries.begin(), tile_boundaries.end());
      }
      return std::make_shared<TiledArray::TiledRange>(ranges.begin(), ranges.end());
    }

  protected:
    std::vector<IndexRange> range_;

    bool initialized_;

//TODO do something for symmetry
//  std::map<std::vector<int>, std::pair<double,bool>> perm_;

    std::tuple<std::string, std::array<size_t, N>, std::array<size_t, N>> index_mapping(const std::string& vars) const {
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
      for (auto i = ss.rbegin(); i != ss.rend(); ++i) {
        if (i != ss.rbegin()) revvars += ",";
        revvars += *i;
      }

      std::array<size_t, N> low, high;
      auto riter = range_.rbegin();
      int n = 0;
      for (auto i = ss.rbegin(); i != ss.rend(); ++i, ++riter, ++n) {
        regex re2("[a-z]+");
        smatch m;
        regex_search(*i, m, re2);
        std::string label = m[0]; // c, x, a, ci, o

        bool find = false;
        low[n] = 0;
        high[n] = 0;
        for (int j = 0; j != riter->nblock(); ++j) {
          const bool matches = riter->range(j).label() == label;
          find |= matches;
          if (!find) low[n] += 1;
          if (!find || matches) high[n] += 1;
        }
      }
      return std::make_tuple(revvars, low, high);
    }

    struct DotProduct {
      using result_type = DataType;
      using first_argument_type = value_type;
      using second_argument_type = value_type;
      result_type operator()() const { return 0.0; }
      const result_type& operator()(const result_type& result) const { return result; }
      void operator()(result_type& result, const result_type& arg) const { result += arg; }
      void operator()(result_type& result, const first_argument_type& first, const second_argument_type& second) const {
        assert(!second.empty());
        // TODO not the best code at the moment
        result += first.conj().dot(second);
      }
    };
    struct Size {
      using result_type = size_t;
      using argument_type = value_type;
      result_type operator()() const { return 0lu; }
      const result_type& operator()(const result_type& result) const { return result; }
      void operator()(result_type& result, const result_type& arg) const { result += arg; }
      void operator()(result_type& result, const argument_type& arg) const { result += arg.size(); }
    };

  public:
    TATensor(const std::vector<IndexRange>& r, const bool initialize = false, /*toggle for symmetry*/const bool dummy = false)
      : BaseArray(madness::World::get_default(), *make_trange(r)), range_(r), initialized_(initialize) {
      static_assert(N != 0, "this should not happen");
      assert(r.size() == N);
      if (initialize)
        this->fill_local(0.0);
    }

    TATensor(const TATensor<DataType,N>& o) : BaseArray(TiledArray::clone_part(o)), range_(o.range_), initialized_(true) {
      static_assert(N != 0, "this should not happen");
    }

    TATensor(TATensor<DataType,N>&& o) : BaseArray(std::move(o)), range_(o.range_), initialized_(true) {
      static_assert(N != 0, "this should not happen");
    }

    virtual void init() { assert(false); }

    void zero() {
      const DataType zero = static_cast<DataType>(0.0);
      for (auto it = begin(); it != end(); ++it)
        if (it->probe())
          get_world().taskq.add([=](value_type& x) { std::fill(x.begin(), x.end(), static_cast<DataType>(0.0)); }, (*it).future());
    }

    void scale(const DataType& a) {
      for (auto it = begin(); it != end(); ++it)
        if (it->probe())
          get_world().taskq.add([=](value_type& x) { x.scale_to(a); }, (*it).future());
    }

    void ax_plus_y(const DataType& a, std::shared_ptr<const TATensor<DataType,N>> o) { ax_plus_y(a, *o); }
    void ax_plus_y(const DataType& a, const TATensor<DataType,N>& o) {
      assert(range_ == o.range_);

      for (auto it = begin(); it != end(); ++it)
        if (it->probe())
          get_world().taskq.add([=](value_type& y, const value_type& x) {
                                  assert(!x.empty());
                                  y.inplace_binary(x, [=](DataType& l, const DataType r) { l += r*a; });
                                }, (*it).future(), o.find(it.ordinal()));
    }

    DataType dot_product(std::shared_ptr<const TATensor<DataType,N>> o) const { return dot_product(*o); }
    DataType dot_product(const TATensor<DataType,N>& o) const {
      assert(range_ == o.range_);

      TiledArray::detail::ReducePairTask<DotProduct> reduce_task(get_world());
      for (auto it = begin(); it != end(); ++it) {
        if (it->probe())
          reduce_task.add((*it).future(), o.find(it.ordinal()));
      }
      return get_world().gop.all_reduce(id(), reduce_task.submit(), DotProduct()).get();
    }

    size_t size_alloc() const {
      TiledArray::detail::ReduceTask<Size> reduce_task(get_world());
      for (auto it = begin(); it != end(); ++it)
        if (it->probe())
          reduce_task.add((*it).future());
      return get_world().gop.all_reduce(id(), reduce_task.submit(), Size()).get();
    }

    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    double rms() const { return std::sqrt(detail::real(dot_product(*this))/size_alloc()); }

    bool initialized() const { return initialized_; }

    void fill_local(const DataType& o) {
      initialized_ = true;
      BaseArray::fill_local(o);
    }

    void init_tile(typename BaseArray::iterator it) {
      madness::Future<value_type> t
        = get_world().taskq.add([](typename BaseArray::range_Type r) {
                                  value_type tile(r);
                                  std::fill(tile.begin(), tile.end(), 0.0);
                                  return tile;
                                }, trange().make_tile_range(it.ordinal()));
      *it = t;
    }

    auto get_local(const std::vector<Index>& index) -> decltype(std::make_pair(true,begin())) {
      assert(index.size() == N);
      // find index and set lo
      std::array<size_t,N> lo_in;
      for (int n = 0; n != N; ++n) {
        auto iter = std::find_if(range_[n].begin(), range_[n].end(), [&](const Index& i) { return i.offset() == index[n].offset(); });
        lo_in[n] = std::accumulate(range_[n].begin(), iter, 0lu, [](size_t n, const Index& i){ return n + i.size(); });
      }
      auto it = begin();
      for ( ; it != end(); ++it) {
        const TiledArray::Range range = trange().make_tile_range(it.ordinal());
        auto lo = range.lobound();
        assert(lo.size() == N);
        if (std::equal(lo_in.rbegin(), lo_in.rend(), lo.begin()))
          break;
      }
      return std::make_pair(it != end(), it);
    }

    std::vector<IndexRange> indexrange() const { return range_; }

    auto operator()(const std::string& vars) const -> decltype(BaseArray::operator()("").block({0},{0})) {
      auto m = index_mapping(vars);
      return BaseArray::operator()(std::get<0>(m)).block(std::get<1>(m), std::get<2>(m));
    }

    auto operator()(const std::string& vars) -> decltype(BaseArray::operator()("").block({0},{0})) {
      auto m = index_mapping(vars);
      return BaseArray::operator()(std::get<0>(m)).block(std::get<1>(m), std::get<2>(m));
    }

    TATensor<DataType,N>& operator=(TATensor<DataType,N>&& o) {
      range_ = o.range_;
      initialized_ = o.initialized_;
      BaseArray::operator=(o);
      return *this;
    }

    TATensor<DataType,N>& operator=(const TATensor<DataType,N>& o) {
      range_ = o.range_;
      initialized_ = o.initialized_;
      BaseArray::operator=(TiledArray::clone_part(o));
      return *this;
    }

    std::shared_ptr<TATensor<DataType,N>> clone() const { return std::make_shared<TATensor<DataType,N>>(range_); }
    std::shared_ptr<TATensor<DataType,N>> copy() const { return std::make_shared<TATensor<DataType,N>>(*this); }

    // TODO temp solution to complex conjugate
    std::shared_ptr<TATensor<DataType,N>> conjg() const {
      std::shared_ptr<TATensor<DataType,N>> out = copy();
      for (auto it = out->begin(); it != out->end(); ++it)
        if (it->probe())
          get_world().taskq.add([=](value_type& x) { blas::conj_n(x.begin(), x.size()); }, (*it).future());
      return out;
    }

    // Dummy function...
    DataType get_scalar() const { assert(false); return 0.0; }
};


template<typename DataType>
class TATensor<DataType,0> : public TiledArray::Array<DataType,1> {
  public:
    using TiledArray::Array<DataType,1>::end;
    using TiledArray::Array<DataType,1>::trange;
    auto begin() -> decltype(TiledArray::Array<DataType,1>::end()) { return end(); }
    auto begin() const -> decltype(TiledArray::Array<DataType,1>::end()) { return end(); }

  protected:
    std::vector<IndexRange> range_;
    DataType data_;

    bool initialized_;

    static std::shared_ptr<TiledArray::TiledRange> make_trange(const std::vector<IndexRange>&) {
      std::vector<TiledArray::TiledRange1> ranges;
      std::vector<size_t> tile_boundaries {0,1};
      ranges.emplace_back(tile_boundaries.begin(), tile_boundaries.end());
      return std::make_shared<TiledArray::TiledRange>(ranges.begin(), ranges.end());
    }

  public:
    TATensor(const std::vector<IndexRange>& r, /*dummy*/bool initialize = true, /*toggle for symmetry*/const bool dummy = false)
      : TiledArray::Array<DataType,1>(madness::World::get_default(), *make_trange(r)), range_(r), initialized_(true) {
    }

    TATensor(const TATensor<DataType,0>& o) : TiledArray::Array<DataType,1>(o), range_(o.range_), data_(o.data_), initialized_(true) {
    }

    TATensor(TATensor<DataType,0>&& o) : TiledArray::Array<DataType,1>(std::move(o)), range_(o.range_), data_(o.data_), initialized_(true) {
    }

    virtual void init() { assert(false); }

    DataType& operator()(const std::string& vars) { assert(vars.empty()); return data_; }
    const DataType& operator()(const std::string& vars) const { assert(vars.empty()); return data_; }

    TATensor<DataType,0>& operator=(TATensor<DataType,0>&& o) {
      range_ = o.range_;
      data_ = o.data_;
      initialized_ = o.initialized_;
      TiledArray::Array<DataType,1>::operator=(o);
      return *this;
    }

    TATensor<DataType,0>& operator=(const TATensor<DataType,0>& o) {
      range_ = o.range_;
      data_ = o.data_;
      initialized_ = o.initialized_;
      TiledArray::Array<DataType,1>::operator=(o);
      return *this;
    }

    std::shared_ptr<TATensor<DataType,0>> clone() const { return std::make_shared<TATensor<DataType,0>>(range_); }
    std::shared_ptr<TATensor<DataType,0>> copy() const { return std::make_shared<TATensor<DataType,0>>(*this); }

    std::vector<IndexRange> indexrange() const { return range_; }

    void zero() { data_ = static_cast<DataType>(0.0); }
    void scale(const DataType& a) { data_ *= a; }

    bool initialized() const { return initialized_; }
    void fill_local(const DataType& i) { data_ = i; }

    DataType get_scalar() const { return data_; }
};

template<typename DataType>
std::ostream& operator<<(std::ostream& os, const TATensor<DataType,0>& o) { os << o(""); return os; }

}}

#endif
