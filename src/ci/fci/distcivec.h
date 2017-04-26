//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: distcivec.h
// Copyright (C) 2011 Toru Shiozaki
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


#ifndef BAGEL_FCI_DISTCIVEC_H
#define BAGEL_FCI_DISTCIVEC_H

#include <bagel_config.h>
#include <src/ci/fci/civec.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/util/parallel/rmawindow.h>

namespace bagel {

template<typename DataType>
class DistCivector : public RMAWindow<DataType> {
  public:
    using DetType = Determinants;
    using LocalizedType = std::false_type;

    using RMAWindow<DataType>::scale;
    using RMAWindow<DataType>::ax_plus_y;
    using RMAWindow<DataType>::dot_product;
    using RMAWindow<DataType>::fence;
    using RMAWindow<DataType>::fence_local;
    using RMAWindow<DataType>::local_data;

  protected:
    mutable std::shared_ptr<const Determinants> det_;

    // global dimension
    size_t lena_;
    size_t lenb_;

    // table for alpha string distribution
    StaticDist dist_;

    // local alpha strings
    size_t astart_;
    size_t aend_;

  public:
    DistCivector(std::shared_ptr<const Determinants> det);
    DistCivector(std::shared_ptr<const Civector<DataType>> civ);

    DistCivector(const DistCivector<DataType>& o) : DistCivector(o.det()) { RMAWindow<DataType>::operator=(o); }
    DistCivector(std::shared_ptr<const DistCivector<DataType>> o) : DistCivector(*o) {}

    // functions required by RMAWindow
    bool is_local(const size_t a) const override;
    std::tuple<size_t, size_t, size_t> locate(const size_t a) const override;
    size_t localsize() const override { return size(); }

    std::shared_ptr<DistCivector<DataType>> clone() const { return std::make_shared<DistCivector<DataType>>(det_); }
    std::shared_ptr<DistCivector<DataType>> copy() const { return std::make_shared<DistCivector<DataType>>(*this); }

    size_t size() const { return lenb_*(aend_-astart_); }
    size_t global_size() const { return lena_*lenb_; }
    size_t lena() const { return lena_; }
    size_t lenb() const { return lenb_; }

    size_t astart() const { return astart_; }
    size_t aend() const { return aend_; }
    size_t asize() const { return aend_ - astart_; }

    const DataType* data() const { return local_data(); }

    void synchronize(const int root = 0) { /* do nothing */ }

    std::shared_ptr<Civector<DataType>> civec() const;
    std::shared_ptr<const Determinants> det() const { return det_; }

    void set_det(std::shared_ptr<const Determinants> d) { det_ = d; }
    void set_local(const size_t la, const size_t lb, const DataType a);

    // utility functions
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    double variance() const { return detail::real(dot_product(*this)) / (lena_*lenb_); }
    double rms() const { return std::sqrt(variance()); }
    void project_out(std::shared_ptr<const DistCivector<DataType>> o) { ax_plus_y(-detail::conj(dot_product(*o)), *o); }

    DataType spin_expectation() const {
      std::shared_ptr<DistCivector<DataType>> S2 = spin();
      return dot_product(*S2);
    }
    std::shared_ptr<DistCivector<DataType>> spin() const;
    void spin_decontaminate(const double thresh = 1.0e-4);

    double orthog(std::list<std::shared_ptr<const DistCivector<DataType>>> c) {
      for (auto& iter : c)
        project_out(iter);
      return normalize();
    }

    double orthog(std::shared_ptr<const DistCivector<DataType>> o) {
      return orthog(std::list<std::shared_ptr<const DistCivector<DataType>>>{o});
    }

    double normalize() {
      const double norm = this->norm();
      const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
      scale(static_cast<DataType>(scal));
      return norm;
    }

    std::shared_ptr<DistCivector<DataType>> transpose() const;

    void print(const double thresh = 0.05) const;
};

template <> void DistCivector<double>::spin_decontaminate(const double);
template <> void DistCivector<std::complex<double>>::spin_decontaminate(const double);

extern template class DistCivector<double>;
extern template class DistCivector<std::complex<double>>;

using DistCivec = DistCivector<double>;
using ZDistCivec = DistCivector<std::complex<double>>;

using DistDvec = Dvector_base<DistCivec>;

}

#endif
