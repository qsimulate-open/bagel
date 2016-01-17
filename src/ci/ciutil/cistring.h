//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cistring.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef BAGEL_CIUTIL_STRINGSPACE_H
#define BAGEL_CIUTIL_STRINGSPACE_H

#include <bitset>
#include <algorithm>
#include <src/util/constants.h>
#include <src/util/parallel/staticdist.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/util/serialization.h>
#include <src/ci/ciutil/cistringmap.h>

namespace bagel {

// Contains all the strings and information for lexical ordering for one particular graph (set of strings)
//   comprised of three subgraphs (one each for RASI, RASII, RASIII)
class CIGraph {
  protected:
    size_t nele_;
    size_t norb_;
    size_t size_;
    std::vector<size_t> weights_;

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & nele_ & norb_ & size_& weights_;
    }

  public:
    CIGraph() { }
    CIGraph(const size_t nele, const size_t norb);

    size_t& weight(const size_t i, const size_t j) { assert(nele_*norb_ > 0); return weights_[i + j*norb_]; }
    const size_t& weight(const size_t i, const size_t j) const { assert(nele_*norb_ > 0); return weights_[i + j*norb_]; }

    size_t size() const { return size_; }

    size_t lexical(const int& start, const int& fence, const std::bitset<nbit__>& abit) const {
      size_t out = 0;

      int k = 0;
      for (int i = start; i < fence; ++i)
        if (abit[i]) out += weight(i-start,k++);
      return out;
    }
};


class CIString_base {
  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int) { }
  public:
    CIString_base() { }
    virtual ~CIString_base() { }
};


template<int N, class Derived>
class CIString_base_impl : public CIString_base {
  protected:
    int norb_;
    int nele_;
    size_t offset_;
    std::vector<std::bitset<nbit__>> strings_;

    std::array<std::pair<int, int>, N> subspace_;
    std::array<std::shared_ptr<CIGraph>, N> graphs_;
    std::shared_ptr<const StaticDist> dist_;

    void init() {
      compute_strings();
    }
    void compute_strings() { static_cast<Derived*>(this)->compute_strings_impl(); }

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template <class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<CIString_base>(*this) << norb_ << nele_ << offset_ << strings_ << subspace_ << graphs_;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<CIString_base>(*this) >> norb_ >> nele_ >> offset_ >> strings_ >> subspace_ >> graphs_;
      const size_t size = std::accumulate(graphs_.begin(), graphs_.end(), 1, [](size_t n, const std::shared_ptr<CIGraph>& i) { return n*i->size(); });
      dist_ = std::make_shared<StaticDist>(size, mpi__->size());
    }

  public:
    CIString_base_impl() : norb_(0), nele_(0), offset_(0) { }
    CIString_base_impl(std::initializer_list<size_t> args) {
      assert(args.size() == 2*N+1);
      auto iter = args.begin();
      for (int i = 0; i != N; ++i) {
        const int a = *iter++;
        const int b = *iter++;
        subspace_[i] = {a, b};
        graphs_[i] = std::make_shared<CIGraph>(a, b);
      }
      // setting to CIString_base
      offset_ = *iter++;
      norb_ = std::accumulate(subspace_.begin(), subspace_.end(), 0, [](int n, const std::pair<int, int>& i) { return n+i.second; });
      nele_ = std::accumulate(subspace_.begin(), subspace_.end(), 0, [](int n, const std::pair<int, int>& i) { return n+i.first; });
      assert(iter == args.end());

      const size_t size = std::accumulate(graphs_.begin(), graphs_.end(), 1, [](size_t n, const std::shared_ptr<CIGraph>& i) { return n*i->size(); });
      dist_ = std::make_shared<StaticDist>(size, mpi__->size());
    }

    // copy construct with an offset update
    CIString_base_impl(const CIString_base_impl<N,Derived>& o, const size_t offset)
      : norb_(o.norb_), nele_(o.nele_), offset_(offset), strings_(o.strings_), subspace_(o.subspace_),
        graphs_(o.graphs_), dist_(o.dist_) { }

    virtual ~CIString_base_impl() { }

    int nele() const { return nele_; }
    int norb() const { return norb_; }
    size_t size() const { return strings_.size(); }
    size_t offset() const { return offset_; }

    bool empty() const { return size() == 0; }

    const std::vector<std::bitset<nbit__>>& strings() const { return strings_; }
    const std::bitset<nbit__>& strings(const size_t i) const { return strings_[i]; }

    std::vector<std::bitset<nbit__>>::iterator begin() { return strings_.begin(); }
    std::vector<std::bitset<nbit__>>::iterator end() { return strings_.end(); }
    std::vector<std::bitset<nbit__>>::const_iterator begin() const { return strings_.cbegin(); }
    std::vector<std::bitset<nbit__>>::const_iterator end() const { return strings_.cend(); }

    size_t size(const int& i) const { return graphs_[i]->size(); }

    std::shared_ptr<const StaticDist> dist() const { return dist_; }

    size_t lexical_zero(const std::bitset<nbit__>& bit) const { return static_cast<const Derived*>(this)->lexical_zero_impl(bit); }
    size_t lexical_offset(const std::bitset<nbit__>& bit) const { return static_cast<const Derived*>(this)->lexical_offset_impl(bit); }

    bool contains(const std::bitset<nbit__>& bit) const { return static_cast<const Derived*>(this)->contains_impl(bit); }
    bool matches(const int i, const int j) const { return static_cast<const Derived*>(this)->matches_impl(i,j); }
    template <typename U>
    bool matches(std::shared_ptr<const U> o) const { return static_cast<const Derived*>(this)->matches_impl(o); }
};


class RASString : public CIString_base_impl<3,RASString> {
  friend class CIString_base_impl<3,RASString>;
  protected:
    void compute_strings_impl();

    bool contains_impl(const std::bitset<nbit__>& bit) const {
      assert(bit.count() == nele_);
      return nholes(bit) == nholes() && nparticles(bit) == nparticles();
    }

    bool matches_impl(const int nh, const int np) const {
      return nh == nholes() && np == nparticles();
    }
    bool matches_impl(const std::shared_ptr<const RASString> o) const {
      return matches_impl(o->nholes(), o->nparticles());
    }

    size_t lexical_offset_impl(const std::bitset<nbit__>& bit) const {
      return lexical_zero(bit) + offset_;
    }

    size_t lexical_zero_impl(const std::bitset<nbit__>& bit) const {
      const size_t r1 = subspace_[0].second;
      const size_t r2 = subspace_[1].second;
      const size_t r3 = subspace_[2].second;

      const size_t n2 = graphs_[1]->size();
      const size_t n1 = graphs_[0]->size();

      size_t out = 0;
      out += graphs_[1]->lexical(r1, r1+r2, bit);
      out += n2 * graphs_[0]->lexical(0, r1, bit);
      out += n2 * n1 * graphs_[2]->lexical(r1+r2, r1+r2+r3, bit);

      return out;
    }

    // helper functions
    int nholes(const std::bitset<nbit__>& bit) const {
      return subspace_[0].second - (bit & (~std::bitset<nbit__>(0ull) >> (nbit__ - subspace_[0].second))).count();
    }
    int nparticles(const std::bitset<nbit__>& bit) const {
      return (bit & (~(~std::bitset<nbit__>(0ull) << subspace_[2].second) << subspace_[0].second + subspace_[1].second)).count();
    }

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<CIString_base_impl<3,RASString>>(*this);
    }

  public:
    RASString() { }
    RASString(const size_t nele1, const size_t norb1, const size_t nele2, const size_t norb2, const size_t nele3, const size_t norb3, const size_t offset = 0);
    RASString(const RASString& o, const size_t offset = 0) : CIString_base_impl<3,RASString>(o, offset) { }

    int nholes() const { return subspace_[0].second - subspace_[0].first; }
    int nele2() const { return nele_ - subspace_[0].first - subspace_[2].first; }
    int nparticles() const { return subspace_[2].first; }

    size_t tag() const { return nholes() + (nparticles() << 8); }

    template <int S> const std::pair<const int, const int> ras() const {
      static_assert(S == 0 || S == 1 || S == 2, "illegal call of RAString::ras");
      return std::get<S>(subspace_);
    }

};


class FCIString : public CIString_base_impl<1,FCIString> {
  friend class CIString_base_impl<1,FCIString>;
  protected:
    void compute_strings_impl();

    bool contains_impl(const std::bitset<nbit__>& bit) const { assert(bit.count() == nele_); return true; }
    bool matches_impl(const int n, const int m) const { return true; }
    bool matches_impl(const std::shared_ptr<const FCIString> o) const { return true; }

    size_t lexical_offset_impl(const std::bitset<nbit__>& bit) const { return lexical(bit)+offset_; }
    size_t lexical_zero_impl(const std::bitset<nbit__>& bit) const { return lexical(bit); }

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<CIString_base_impl<1,FCIString>>(*this);
    }

  public:
    FCIString() { }
    FCIString(const size_t nele1, const size_t norb1, const size_t offset = 0);
    FCIString(const FCIString& o, const size_t offset = 0) : CIString_base_impl<1,FCIString>(o, offset) { }

    size_t lexical(const std::bitset<nbit__>& bit) const {
      assert(contains(bit));
      return graphs_[0]->lexical(0, norb_, bit);
    }
};

using FCIString_base = CIString_base_impl<1,FCIString>;
using RASString_base = CIString_base_impl<3,RASString>;

}

extern template class bagel::CIString_base_impl<1,bagel::FCIString>;
extern template class bagel::CIString_base_impl<3,bagel::RASString>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::CIGraph)
BOOST_CLASS_EXPORT_KEY(bagel::CIString_base)
BOOST_CLASS_EXPORT_KEY(bagel::RASString_base)
BOOST_CLASS_EXPORT_KEY(bagel::FCIString_base)
BOOST_CLASS_EXPORT_KEY(bagel::RASString)
BOOST_CLASS_EXPORT_KEY(bagel::FCIString)

#endif
