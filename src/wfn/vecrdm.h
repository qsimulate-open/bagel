//
// BAGEL - Parallel electron correlation program.
// Filename: vecrdm.h
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

#ifndef __SRC_WFN_VECRDM_H
#define __SRC_WFN_VECRDM_H

#include <src/wfn/rdm.h>

namespace bagel {

template<int N>
class VecRDM {
  protected:
    typename std::map<std::pair<int,int>, std::shared_ptr<RDM<N>>> data_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & data_; }

  public:
    VecRDM() { }
    VecRDM(const VecRDM& o) {
      for (auto& i : o.data_)
        data_.emplace(i.first, i.second->copy());
    }

    size_t size() const { return data_.size(); }

    typename std::map<std::pair<int,int>, std::shared_ptr<RDM<N>>>::iterator begin() { return data_.begin(); }
    typename std::map<std::pair<int,int>, std::shared_ptr<RDM<N>>>::iterator end() { return data_.end(); }
    typename std::map<std::pair<int,int>, std::shared_ptr<RDM<N>>>::const_iterator begin() const { return data_.cbegin(); }
    typename std::map<std::pair<int,int>, std::shared_ptr<RDM<N>>>::const_iterator end() const { return data_.cend(); }
    typename std::map<std::pair<int,int>, std::shared_ptr<RDM<N>>>::const_iterator cbegin() const { return data_.cbegin(); }
    typename std::map<std::pair<int,int>, std::shared_ptr<RDM<N>>>::const_iterator cend() const { return data_.cend(); }

    // adding elements (i and j are bra and ket state indices)
    void emplace(const int i, const int j, std::shared_ptr<RDM<N>> d) {
      // if this key is present, the element is deleted
      auto key = std::make_pair(i, j);
      if (data_.find(key) != data_.end())
        data_.erase(key);
      data_.emplace(key, d);
    }
    // special function for diagonal RDM
    void emplace(const int i, std::shared_ptr<RDM<N>> d) { emplace(i, i, d); }

    // get RDMs
    std::shared_ptr<RDM<N>> at(const int i) { return at(i, i); }
    std::shared_ptr<RDM<N>> at(const int i, const int j) { return data_.at(std::make_pair(i, j)); }
    std::shared_ptr<const RDM<N>> at(const int i) const { return at(i, i); }
    std::shared_ptr<const RDM<N>> at(const int i, const int j) const { return data_.at(std::make_pair(i, j)); }

};

}

#endif
