//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cistringmap.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __SRC_CIUTIL_CISTRINGMAP_H
#define __SRC_CIUTIL_CISTRINGMAP_H

namespace bagel {

struct DetMap {
  public:
    size_t target;
    int sign;
    size_t source;
    int ij;

    DetMap() { }
    DetMap(size_t t, int si, size_t s, int o) : target(t), sign(si), source(s), ij(o) {}

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & target & source & sign; }
};


class StringMap {
  protected:
    std::vector<std::vector<DetMap>> data_;
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & data_; }
  public:
    StringMap() { }
    StringMap(const int norb) : data_(norb) { }
    StringMap(const StringMap& o) : data_(o.data_) { }

    void clear() { data_.clear(); }
    void resize(const size_t n) { data_.resize(n); }

    void reserve(const size_t n) {
      for (auto& i : data_)
        i.reserve(n);
    }

    void shrink_to_fit() {
      for (auto& i : data_)
        i.shrink_to_fit();
    }

    std::vector<DetMap>& operator[](const size_t i) { return data_[i]; }
    const std::vector<DetMap>& operator[](const size_t i) const { return data_[i]; }

    const std::vector<DetMap>& data(const size_t i) const { assert(i < data_.size()); return data_[i]; }

    void insert(const std::vector<std::vector<DetMap>>& inp) {
      assert(data_.size() == inp.size());
      auto j = inp.begin();
      for (auto& i : data_) {
        i.insert(i.end(), j->begin(), j->end());
        ++j;
      }
    }
    void insert(const std::shared_ptr<const StringMap>& o) { insert(o->data_); }

    size_t size() const { return std::accumulate(data_.begin(), data_.end(), 0, [](size_t n, const std::vector<DetMap>& o) { return n+o.size(); }); }

    std::shared_ptr<StringMap> get_minus() const {
      auto out = std::make_shared<StringMap>(*this);
      for (auto& i : out->data_)
        for (auto& j : i)
          j.sign = -j.sign;
      return out;
    }
};


template <class StringType>
class DetMapBlock_base {
  protected:
    size_t offset_;
    std::shared_ptr<const StringType> source_space_;
    std::vector<DetMap> phis_;
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & offset_ & source_space_ & phis_; }

  public:
    DetMapBlock_base() { }
    DetMapBlock_base(const size_t o, std::shared_ptr<const StringType> sp, std::vector<DetMap>&& p) : offset_(o), source_space_(sp), phis_(std::move(p)) {}


    std::vector<DetMap>::const_iterator begin() const { return phis_.begin(); }
    std::vector<DetMap>::const_iterator end() const { return phis_.end(); }
    size_t size() const { return phis_.size(); }
    size_t offset() const { return offset_; }
    std::shared_ptr<const StringType> source_space() const { return source_space_; }
};

}
#endif
