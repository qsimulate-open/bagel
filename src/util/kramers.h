//
// BAGEL - Parallel electron correlation program.
// Filename: kramers.h
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

#ifndef __SRC_UTIL_KRAMERS_H
#define __SRC_UTIL_KRAMERS_H

#include <map>
#include <memory>
#include <bitset>
#include <string>
#include <sstream>
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

    std::bitset<N> tag() const { return tag_; }

    KTag<N> perm(const std::array<int,N>& o) const {
      std::bitset<N> out;
      for (int i = 0; i != N; ++i)
        out[i] = tag_[o[i]];
      return KTag<N>(out);
    }

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
    std::map<std::array<int,N>, double> perm_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      ar & data_ & perm_;
    }

  public:
    Kramers() { }
    void emplace(const KTag<N>& t, std::shared_ptr<Type> o) { assert(!exist(t)); data_.emplace(t, o); }

    size_t size() const { return data_.size(); }

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

    void emplace_perm(const std::array<int,N>& o, double a) { perm_.emplace(o,a); }

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
              sort_indices(/*sort info*/i.first, /*fac*/i.second, /*fac2*/0.0, j.second->data(), out->data(), dim);
              break;
            }
          }
          if (found) break;
        }
        return out;
      }
    }
};

}

#endif
