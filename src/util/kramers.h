//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: kramers.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_UTIL_KRAMERS_H
#define __SRC_UTIL_KRAMERS_H

#include <map>
#include <memory>
#include <bitset>
#include <string>
#include <sstream>
#include <array>
#include <initializer_list>
#include <cassert>
#include <src/util/serialization.h>
#include <src/util/prim_op_var.h>

namespace bagel {

template<int N>
class KTag {
  protected:
    std::bitset<N> tag_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      ar & tag_;
    }

  public:
    KTag() {}
    KTag(const std::bitset<N>& o) : tag_(o) { }
    KTag(const std::string& o) : tag_(o) { }
    KTag(const unsigned long o) : tag_(o) { }
    KTag(std::initializer_list<int> il) {
      std::stringstream ss;
      for (auto& i : il) {
        assert(i/2 == 0);
        ss << i;
      }
      tag_ = std::bitset<N>(ss.str());
    }
    KTag(const std::vector<bool>& o) {
      assert(o.size() == N);
      for (int i = 0; i != N; ++i)
        tag_[N-1-i] = o[i];
    }
    KTag(const std::array<bool,N>& o) {
      for (int i = 0; i != N; ++i)
        tag_[N-1-i] = o[i];
    }

    std::bitset<N> tag() const { return tag_; }

    std::vector<bool> vec() const {
      std::vector<bool> o(N);
      for (int i = 0; i != N; ++i)
        o[i] = tag_[N-1-i];
      return o;
    }

    KTag<N> perm(const std::vector<int>& o) const {
      assert(o.size() == N);
      std::bitset<N> out;
      for (int i = 0; i != N; ++i)
        out[N-i-1] = tag_[N-o[i]-1];
      return KTag<N>(out);
    }

    KTag<N>& operator=(const KTag<N>& o) { tag_ = o.tag_; return *this; }
    bool operator==(const KTag<N>& o) const { return tag_.to_string() == o.tag_.to_string(); }
    bool operator!=(const KTag<N>& o) const { return !(*this == o); }
    bool operator<(const KTag<N>& o) const { return tag_.to_string() < o.tag_.to_string(); }
    bool operator<=(const KTag<N>& o) const { return *this < o || *this == o; }
    bool operator>(const KTag<N>& o) const { return !(*this <= o); }
    bool operator>=(const KTag<N>& o) const { return !(*this < o); }
};

template<int N, int M>
KTag<N+M> merge(const KTag<N>& a, const KTag<M>& b) { return KTag<N+M>(a.tag().to_string()+b.tag().to_string()); }

template<int N, class Type>
class Kramers {
  protected:
    std::map<KTag<N>, std::shared_ptr<Type>> data_;
    std::map<std::vector<int>, std::pair<double,bool>> perm_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      ar & data_ & perm_;
    }

  public:
    Kramers() { }
    Kramers(const Kramers<N,Type>& o) : perm_(o.perm_) {
      for (auto& i : o.data_)
        data_.emplace(i.first, i.second->copy());
    }
    void emplace(const KTag<N>& t, std::shared_ptr<Type> o) { assert(!exist(t)); data_.emplace(t, o); }

    void zero() {
      for (auto& i : data_)
        i.second->zero();
    }

    std::shared_ptr<Kramers<N,Type>> copy() const { return std::make_shared<Kramers<N,Type>>(*this); }

    std::list<std::vector<bool>> listkramers() const {
      std::list<std::vector<bool>> out;
      for (auto& i : data_) {
        std::vector<bool> v = i.first.vec();
        if (std::find(out.begin(), out.end(), v) == out.end())
          out.push_back(v);
      }
      return out;
    }

    size_t size() const { return data_.size(); }
    const std::map<std::vector<int>, std::pair<double,bool>>& perm() const { return perm_; }

    std::shared_ptr<Type>& at(const KTag<N>& i) { return data_.at(i); }
    std::shared_ptr<const Type> at(const KTag<N>& i) const { return data_.at(i); }

    bool exist(const KTag<N>& t) const { return data_.find(t) != data_.end(); }

    std::shared_ptr<Type>& operator[](const KTag<N>& t) { return data_[t]; }
    typename std::map<KTag<N>, std::shared_ptr<Type>>::iterator find(const KTag<N>& t) { return data_.find(t); }

    typename std::map<KTag<N>, std::shared_ptr<Type>>::iterator begin() { return data_.begin(); }
    typename std::map<KTag<N>, std::shared_ptr<Type>>::iterator end()   { return data_.end(); }
    typename std::map<KTag<N>, std::shared_ptr<Type>>::const_iterator begin() const { return data_.cbegin(); }
    typename std::map<KTag<N>, std::shared_ptr<Type>>::const_iterator end()   const { return data_.cend(); }
    typename std::map<KTag<N>, std::shared_ptr<Type>>::const_iterator cbegin() const { return data_.cbegin(); }
    typename std::map<KTag<N>, std::shared_ptr<Type>>::const_iterator cend()   const { return data_.cend(); }

    void add(const KTag<N>& t, std::shared_ptr<Type> o) {
      if (exist(t))
        *at(t) += *o;
      else
        emplace(t, o);
    }
    void add(const KTag<N>& t, std::shared_ptr<const Type> o) {
      if (exist(t))
        *at(t) += *o;
      else
        emplace(t, o->copy());
    }

    void emplace_perm(const std::vector<int>& o, double a, bool b = false) { assert(N == o.size()); perm_.emplace(o, std::make_pair(a,b)); }
    void set_perm(const std::map<std::vector<int>, std::pair<double,bool>>& o) { perm_ = o; }

    // find the right permutation and sort indices
    std::shared_ptr<const Type> get_data(const KTag<N>& tag) const {
      if (exist(tag)) {
        return at(tag);
      } else {
        std::shared_ptr<Type> out;
        for (auto& i : perm_) {
          bool found = false;
          for (auto& j : data_) {
            if (tag == j.first.perm(i.first)) {
              found = true;
              out = j.second->clone();
              std::array<int,N> dim;
              for (auto& d : dim)
                d = std::lround(std::pow(out->size(), 1.0/N));
              assert(std::pow(dim[0], N) == out->size());
              std::array<int,N> p;
              std::copy(i.first.begin(), i.first.end(), p.data());
              sort_indices(/*sort info*/p, /*fac*/i.second.first, /*fac2*/0.0, j.second->data(), out->data(), dim);
              break;
            }
          }
          if (found) {
            if (i.second.second)
              blas::conj_n(out->data(), out->size());
            break;
          }
        }
        return out;
      }
    }

    void print() const {
      for (auto& i : data_) {
        std::cout << " tag: " << i.first.tag() << std::endl;
        i.second->print();
      }
    }
};

}

#endif
