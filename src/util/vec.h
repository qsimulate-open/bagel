//
// BAGEL - Parallel electron correlation program.
// Filename: vec.h
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

#ifndef __SRC_UTIL_VEC_H
#define __SRC_UTIL_VEC_H

#include <map>
#include <src/util/serialization.h>

namespace bagel {

// Type has to have a copy() function

template<class Type>
class Vec {
  protected:
    typename std::map<std::pair<int,int>, std::shared_ptr<Type>> data_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & data_; }

  public:
    Vec() { }
    Vec(const Vec<Type>& o) {
      for (auto& i : o.data_)
        data_.emplace(i.first, i.second->copy());
    }

    size_t size() const { return data_.size(); }

    typename std::map<std::pair<int,int>, std::shared_ptr<Type>>::iterator begin() { return data_.begin(); }
    typename std::map<std::pair<int,int>, std::shared_ptr<Type>>::iterator end() { return data_.end(); }
    typename std::map<std::pair<int,int>, std::shared_ptr<Type>>::const_iterator begin() const { return data_.cbegin(); }
    typename std::map<std::pair<int,int>, std::shared_ptr<Type>>::const_iterator end() const { return data_.cend(); }
    typename std::map<std::pair<int,int>, std::shared_ptr<Type>>::const_iterator cbegin() const { return data_.cbegin(); }
    typename std::map<std::pair<int,int>, std::shared_ptr<Type>>::const_iterator cend() const { return data_.cend(); }

    // adding elements (i and j are bra and ket state indices)
    void emplace(const int i, const int j, std::shared_ptr<Type> d) {
      // if this key is present, the element is deleted
      auto key = std::make_pair(i, j);
      if (data_.find(key) != data_.end())
        data_.erase(key);
      data_.emplace(key, d);
    }
    // special function for diagonal RDM
    void emplace(const int i, std::shared_ptr<Type> d) { emplace(i, i, d); }

    // get RDMs
    std::shared_ptr<Type> at(const int i) { return at(i, i); }
    std::shared_ptr<Type> at(const int i, const int j) { return data_.at(std::make_pair(i, j)); }
    std::shared_ptr<const Type> at(const int i) const { return at(i, i); }
    std::shared_ptr<const Type> at(const int i, const int j) const { return data_.at(std::make_pair(i, j)); }

    bool exist(const int i, const int j) const { return data_.find(std::make_pair(i, j)) != data_.end(); }
};

}

#endif
