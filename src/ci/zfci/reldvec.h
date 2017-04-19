//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reldvec.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_ZFCI_RELDVEC_H
#define __SRC_ZFCI_RELDVEC_H

#include <src/ci/fci/dvec.h>
#include <src/ci/zfci/relspace.h>

namespace bagel {

template <typename DataType>
class RelDvector {
  protected:
    using MapType = std::map<std::pair<int,int>, std::shared_ptr<Dvector<DataType>>>;
    using EleType = std::pair<std::pair<int,int>, std::shared_ptr<Dvector<DataType>>>;

    MapType dvecs_;
    std::shared_ptr<const Space_base> space_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & dvecs_ & space_;
    }

  public:
    RelDvector() { }
    // make an empty Dvec
    RelDvector(std::shared_ptr<const Space_base> space, const size_t ij);
    RelDvector(const MapType& o, std::shared_ptr<const Space_base> space) : dvecs_(o), space_(space) { }
    RelDvector(const RelDvector<DataType>& o);
    RelDvector(RelDvector<DataType>&& o) : dvecs_(o.dvecs_), space_(o.space_) { }
    // combines (opposite of split())
    RelDvector(const std::vector<std::shared_ptr<RelDvector<DataType>>>& o);

    std::shared_ptr<RelDvector<DataType>> clone() const;
    std::shared_ptr<RelDvector<DataType>> copy() const;

    std::shared_ptr<RelDvector<DataType>> extract_state(const std::vector<int> input) const;

    std::shared_ptr<Dvector<DataType>> find(int a, int b) { return dvecs_.at({a, b}); }
    std::shared_ptr<const Dvector<DataType>> find(int a, int b) const { return dvecs_.at({a, b}); }

    std::shared_ptr<const Space_base> space() const { return space_; }
    MapType dvecs() { return dvecs_; }
    const MapType dvecs() const { return dvecs_; }

    void set_data(const int istate, std::shared_ptr<const RelDvector<DataType>> o);

    void zero();

    size_t size() const;
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    double variance() const { return detail::real(dot_product(*this)) / size(); }
    double rms() const { return std::sqrt(variance()); }

    DataType dot_product(std::shared_ptr<const RelDvector<DataType>> o) const { return dot_product(*o); }
    DataType dot_product(const RelDvector<DataType>& o) const;

    void ax_plus_y(const DataType a, std::shared_ptr<const RelDvector<DataType>> o) { ax_plus_y(a, *o); }
    void ax_plus_y(const DataType a, const RelDvector<DataType>& o);

    std::vector<std::shared_ptr<const RelDvector<DataType>>> split() const { return split(0, dvecs_.begin()->second->ij()); }
    std::vector<std::shared_ptr<const RelDvector<DataType>>> split(const int nstart, const int nend) const;

    std::vector<std::shared_ptr<const RelDvector<DataType>>> dvec(const std::vector<int>& conv) const;

    void project_out(std::shared_ptr<const RelDvector<DataType>> o) { ax_plus_y(-dot_product(*o), *o); }
    void scale(const DataType& a);

    double orthog(std::list<std::shared_ptr<const RelDvector<DataType>>> c);
    double orthog(std::shared_ptr<const RelDvector<DataType>> o);

    double normalize();
    void print(double thresh) const;

    void synchronize();
};

using RelDvec = RelDvector<double>;
using RelZDvec = RelDvector<std::complex<double>>;

}

extern template class bagel::RelDvector<double>;
extern template class bagel::RelDvector<std::complex<double>>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RelDvector<double>)
BOOST_CLASS_EXPORT_KEY(bagel::RelDvector<std::complex<double>>)

#endif
