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
  private:
    TiledArray::Array<DataType, N, LazyTensor>* owner_;
    std::array<size_t, N> index_;
    Gen* gen_;

  public:
    using value_type = DataType;
    using eval_type = TiledArray::Tensor<DataType>;
    using range_type = typename eval_type::range_type;

    template <typename Archive>
    void serialize(Archive& ar) {
      assert(false);
    }

  public:
    LazyTensor(TiledArray::Array<DataType, N, LazyTensor>* o, const std::array<size_t, N>& i, Gen* gen)
     : owner_(o), index_(i), gen_(gen) {
    }

    LazyTensor(const LazyTensor& other) : owner_(other.owner_), index_(other.index_), gen_(other.gen_) { }
    LazyTensor() { }

    LazyTensor& operator=(const LazyTensor& other) {
      owner_ = other.owner_;
      index_ = other.index_;
      gen_   = other.gen_;
      return *this;
    }

    operator eval_type () const {
      eval_type tile(owner_->trange().make_tile_range(index_));
      auto* ptr = tile.data();
      for (auto& i : tile.range())
        *ptr++ = (*gen_)(i);
      return tile;
    }
};


// assuming that the diagonal element is real (which is usually right)
// TODO CAUTION -- column-row major convention.
class Diag4Gen {
  protected:
    std::vector<double> eig_;
    const int o0_, o1_, o2_, o3_;
  public:
    // 0 and 1 are virtual orbitals; 2 and 3 are closed orbitals
    Diag4Gen(const std::vector<double>& eig, const int off0, const int off1, const int off2, const int off3)
     : eig_(eig), o0_(off0), o1_(off1), o2_(off2), o3_(off3) {
    }
    template <typename Index>
    double operator()(const Index& i) {
      // TODO index reversed due to column-row major convention
      return eig_[o0_+i[3]] + eig_[o1_+i[2]] - eig_[o2_+i[1]] - eig_[o3_+i[0]];
    }
};


// TODO CAUTION -- column-row major convention.
class DenomAACC {
  protected:
    std::vector<double> eig_;
    const int o0_, o1_, o2_, o3_;
  public:
    // 0 and 1 are virtual orbitals; 2 and 3 are closed orbitals
    DenomAACC(const std::vector<double>& eig, const int off0, const int off1, const int off2, const int off3)
     : eig_(eig), o0_(off0), o1_(off1), o2_(off2), o3_(off3) {
    }
    template <typename Index>
    double operator()(const Index& i) {
      // TODO index reversed due to column-row major convention
      return 1.0 / (eig_[o0_+i[3]] + eig_[o1_+i[2]] - eig_[o2_+i[1]] - eig_[o3_+i[0]]);
    }
};


template<typename DataType, int N, typename Gen>
class LazyTATensor : public TiledArray::Array<DataType,N,LazyTensor<DataType,N,Gen>> {
  protected:
    using Lazy = LazyTensor<DataType,N,Gen>;
    Gen gen_;

  public:
    LazyTATensor(const std::vector<IndexRange>& r, Gen gen)
      : TiledArray::Array<DataType,N,Lazy>(madness::World::get_default(), *TATensor<DataType,N>::make_trange(r)), gen_(gen) {

      for (auto& t : this->trange().tiles())
        if (this->is_local(t)) {
          std::array<size_t, N> index;
          std::copy(t.begin(), t.end(), index.begin());
          madness::Future<typename LazyTATensor<DataType,N,Gen>::value_type> tile(Lazy(this, index, &gen_));
          this->set(t, tile);
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
      for (auto i = ss.rbegin(); i != ss.rend(); ++i) {
        if (i != ss.rbegin()) revvars += ",";
        revvars += *i;
      }
      return TiledArray::Array<DataType,N,Lazy>::operator()(revvars);
    }
};


}}

#endif
