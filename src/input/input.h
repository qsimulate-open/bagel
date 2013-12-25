//
// BAGEL - Parallel electron correlation program.
// Filename: input.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_INPUT_INPUT_H
#define __SRC_INPUT_INPUT_H

#include <boost/property_tree/ptree.hpp>
#include <array>
#include <sstream>
#include <memory>
#include <vector>
#include <src/util/string.h>

namespace bagel {

// rename as you like
class PTree;
class PTreeIterator {
  private:
    boost::property_tree::ptree::const_iterator current;
  public:
    PTreeIterator() {}
    explicit PTreeIterator(boost::property_tree::ptree::const_iterator curr): current(curr) {}

    //Dereference operator - return the current node's data.
    const std::shared_ptr<const PTree> operator*();

    //Prefix returns by reference.
    PTreeIterator& operator++() { ++current; return *this; }
    PTreeIterator& operator--() { --current; return *this; }

    //Postfix should be implemented in terms of prefix operators
    PTreeIterator operator++(int) { PTreeIterator out = *this; ++*this; return out; }
    PTreeIterator operator--(int) { PTreeIterator out = *this; --*this; return out; }

    bool operator==(const PTreeIterator& o) const { return current == o.current; }
    bool operator!=(const PTreeIterator& o) const { return current != o.current; }

};
class PTreeReverseIterator {
  private:
    boost::property_tree::ptree::const_reverse_iterator current;
  public:
    PTreeReverseIterator() {}
    explicit PTreeReverseIterator(boost::property_tree::ptree::const_reverse_iterator curr): current(curr) {}

    //Dereference operator - return the current node's data.
    const std::shared_ptr<const PTree> operator*();

    //Prefix returns by reference.
    PTreeReverseIterator& operator++() { ++current; return *this; }
    PTreeReverseIterator& operator--() { --current; return *this; }

    //Postfix should be implemented in terms of prefix operators
    PTreeReverseIterator operator++(int) { PTreeReverseIterator out = *this; ++*this; return out; }
    PTreeReverseIterator operator--(int) { PTreeReverseIterator out = *this; --*this; return out; }

    bool operator==(const PTreeReverseIterator& o) const { return current == o.current; }
    bool operator!=(const PTreeReverseIterator& o) const { return current != o.current; }

};

class PTree {
  protected:
    boost::property_tree::ptree data_;
    std::string key_;

  public:
    PTree() : data_() {}

    PTree(const boost::property_tree::ptree& i, const std::string key) : data_(i), key_(key) { }

    PTree(const PTree& o) : data_(o.data_), key_(o.key_) { }

    PTree(const std::string& input);

    std::shared_ptr<PTree> get_child(const std::string& key) const {
      return std::make_shared<PTree>(data_.get_child(key), key);
    }

    std::shared_ptr<PTree> get_child_optional(const std::string& key) const {
      auto out = data_.get_child_optional(key);
      return out ? std::make_shared<PTree>(*out, key) : std::shared_ptr<PTree>();
    }

    template<typename T> T get(const std::string s) const { return data_.get<T>(s); }
    template<typename T> T get(const std::string s, const T& t) const { return data_.get<T>(s, t); }

    void add_child(const std::string s, std::shared_ptr<PTree> ch) { data_.add_child(s, ch->data_); }
    template<typename T> void put(const std::string s, const T& o) { data_.put<T>(s, o); }
    template<typename T> void push_back(const T& o) {
      assert(typeid(T) != typeid(std::shared_ptr<PTree>)); // there is a specialization for shared_ptr<PTree>
      boost::property_tree::ptree ch;
      ch.put("", lexical_cast<std::string>(o));
      data_.push_back(std::make_pair("", ch));
    }

    template<typename T> std::vector<T> get_vector(const std::string s, const int nexpected = 0) const;
    template<typename T, int N> std::array<T,N> get_array(const std::string s) const;

    void erase(const std::string key) { data_.erase(key); }

    std::string data() const { return data_.data(); }
    std::string key() const { return key_; }

    size_t size() const;

    PTreeIterator begin() const;
    PTreeIterator end()   const;
    PTreeReverseIterator rbegin() const;
    PTreeReverseIterator rend()   const;

    void print() const;

    // static function to read basis files
    static std::shared_ptr<const PTree> read_basis(std::string name);
};

template <> void PTree::push_back<std::shared_ptr<PTree>>(const std::shared_ptr<PTree>& pt);

template<typename T> std::vector<T> PTree::get_vector(const std::string key, const int nexpected) const {
  std::vector<T> out;
  auto tmp = get_child(key);
  if ( (nexpected > 0) && (tmp->size() != nexpected) ) {
    std::stringstream err;
    err << "Unexpected number of elements in vector " << key << ". Expected: " << nexpected << ", received: " << tmp->size();
    throw std::runtime_error(err.str());
  }
  for (auto& i : *tmp)
    out.push_back(lexical_cast<T>(i->data()));
  return out;
}

template<typename T, int N> std::array<T,N> PTree::get_array(const std::string key) const {
  std::array<T,N> out;
  auto tmp = get_child(key);
  if (tmp->size() != N) {
    std::stringstream err;
    err << "Unexpected number of elements in array " << key << ". Expected: " << N << ", received: " << tmp->size();
    throw std::runtime_error(err.str());
  }
  int n = 0;
  for (auto& i : *tmp)
    out[n++] = lexical_cast<T>(i->data());
  return out;
}


}

#endif
