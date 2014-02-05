//
// BAGEL - Parallel electron correlation program.
// Filename: ciutil/stringspace.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef BAGEL_CIUTIL_STRINGSPACE_H
#define BAGEL_CIUTIL_STRINGSPACE_H

#include <vector>
#include <array>
#include <memory>
#include <bitset>
#include <algorithm>

#include <src/util/constants.h>
#include <src/parallel/staticdist.h>
#include <src/parallel/mpi_interface.h>
#include <src/util/serialization.h>

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

    const size_t size() const { return size_; }

    const size_t lexical(const int& start, const int& fence, const std::bitset<nbit__>& abit) const {
      size_t out = 0;

      int k = 0;
      for (int i = start; i < fence; ++i)
        if (abit[i]) out += weight(i-start,k++);
      return out;
    }
};


// dummy base class for serialization
namespace detail {
class BaseClass {
  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive&, const unsigned int) { }
  public:
    BaseClass() { }
    virtual ~BaseClass() { }
};
}

template<int N>
class StringSpace_base : public detail::BaseClass {
  protected:
    std::array<std::pair<int, int>, N> subspace_;

    std::vector<std::bitset<nbit__>> strings_;

    std::array<std::shared_ptr<CIGraph>, N> graphs_;

    std::shared_ptr<const StaticDist> dist_;

    int norb_;
    int nele_;

    size_t offset_;

    void init() {
      compute_strings();
    }
    virtual void compute_strings() = 0;

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template <class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<detail::BaseClass>(*this) << subspace_ << strings_ << graphs_;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<detail::BaseClass>(*this) >> subspace_ >> strings_ >> graphs_;
      const size_t size = std::accumulate(graphs_.begin(), graphs_.end(), 1, [](size_t n, const std::shared_ptr<CIGraph>& i) { return n*i->size(); });
      dist_ = std::make_shared<StaticDist>(size, mpi__->size());
    }

  public:
    StringSpace_base() { }
    StringSpace_base(std::initializer_list<size_t> args) {
      assert(args.size() == 2*N+1);
      auto iter = args.begin();
      for (int i = 0; i != N; ++i) {
        const int a = *iter++;
        const int b = *iter++;
        subspace_[i] = std::make_pair(a, b);
        graphs_[i] = std::make_shared<CIGraph>(a, b);
      }
      offset_ = *iter++;
      assert(iter == args.end());

      const size_t size = std::accumulate(graphs_.begin(), graphs_.end(), 1, [](size_t n, const std::shared_ptr<CIGraph>& i) { return n*i->size(); });
      dist_ = std::make_shared<StaticDist>(size, mpi__->size());
      norb_ = std::accumulate(subspace_.begin(), subspace_.end(), 0, [](int n, const std::pair<int, int>& i) { return n+i.second; });
      nele_ = std::accumulate(subspace_.begin(), subspace_.end(), 0, [](int n, const std::pair<int, int>& i) { return n+i.first; });
    }

    virtual ~StringSpace_base() { }

    int nele() const { return nele_; }
    int norb() const { return norb_; }

    template <int subspace> const std::pair<const int, const int> ras() const { return std::get<subspace>(subspace_); }

    size_t size() const { return strings_.size(); }
    size_t offset() const { return offset_; }

    size_t size(const int& i) const { return graphs_[i]->size(); }

    const std::vector<std::bitset<nbit__>>& strings() const { return strings_; }
    const std::bitset<nbit__>& strings(const size_t i) const { return strings_[i]; }

    std::vector<std::bitset<nbit__>>::iterator begin() { return strings_.begin(); }
    std::vector<std::bitset<nbit__>>::iterator end() { return strings_.end(); }
    std::vector<std::bitset<nbit__>>::const_iterator begin() const { return strings_.cbegin(); }
    std::vector<std::bitset<nbit__>>::const_iterator end() const { return strings_.cend(); }

    std::shared_ptr<const StaticDist> dist() const { return dist_; }

    virtual size_t lexical_zero(const std::bitset<nbit__>& bit) const = 0;
    virtual size_t lexical_offset(const std::bitset<nbit__>& bit) const = 0;
};


class StringSpace : public StringSpace_base<3> {
  protected:
    void compute_strings() override;

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<StringSpace_base>(*this);
    }

  public:
    StringSpace() { }
    StringSpace(const size_t nele1, const size_t norb1, const size_t nele2, const size_t norb2, const size_t nele3, const size_t norb3, const size_t offset = 0);

    int nholes() const { return subspace_[0].second - subspace_[0].first; }
    int nele2() const { return nele_ - subspace_[0].first - subspace_[2].first; }
    int nparticles() const { return subspace_[2].first; }

    size_t lexical_offset(const std::bitset<nbit__>& bit) const override {
      return lexical_zero(bit) + offset_;
    }

    size_t lexical_zero(const std::bitset<nbit__>& bit) const override {
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
};


class FCIString : public StringSpace_base<1> {
  protected:
    void compute_strings() override;

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<StringSpace_base>(*this);
    }

  public:
    FCIString() { }
    FCIString(const size_t nele1, const size_t norb1, const size_t offset = 0);

    size_t lexical(const std::bitset<nbit__>& bit) const {
      return graphs_[0]->lexical(0, norb_, bit);
    }
    size_t lexical_offset(const std::bitset<nbit__>& bit) const override { return lexical(bit)+offset_; }
    size_t lexical_zero(const std::bitset<nbit__>& bit) const override { return lexical(bit); }

};

}

extern template class bagel::StringSpace_base<1>;
extern template class bagel::StringSpace_base<3>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::CIGraph)
BOOST_CLASS_EXPORT_KEY(bagel::detail::BaseClass)
BOOST_CLASS_EXPORT_KEY(bagel::StringSpace_base<1>)
BOOST_CLASS_EXPORT_KEY(bagel::StringSpace_base<3>)
BOOST_CLASS_EXPORT_KEY(bagel::StringSpace)
BOOST_CLASS_EXPORT_KEY(bagel::FCIString)

namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<detail::BaseClass, T>::value>::type> {
    typedef detail::BaseClass type;
  };
}

#endif
