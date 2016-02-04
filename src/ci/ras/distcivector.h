//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/civector.h
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


#ifndef __BAGEL_RAS_RASDISTCIVECTOR_H
#define __BAGEL_RAS_RASDISTCIVECTOR_H

#include <list>
#include <src/ci/ras/civector_base.h>
#include <src/ci/ras/apply_block.h>
#include <src/ci/fci/dvector_base.h>

namespace bagel {


template <typename DataType> class RASCivector;
template <typename DataType> using DistCIBlock = DistCIBlock_alloc<DataType, RASString>;

template <typename DataType>
class DistRASCivector : public RASCivector_base<DistCIBlock<DataType>> {
  public: using DetType = RASDeterminants;
  public: using RBlock = DistCIBlock<DataType>;
  public: using LocalizedType = std::false_type;

  protected:
    using RASCivector_base<DistCIBlock<DataType>>::blocks_;
    using RASCivector_base<DistCIBlock<DataType>>::det_;

    // global array tag
    int ga_;

    size_t global_size_;

  public:
    DistRASCivector(std::shared_ptr<const RASDeterminants> det);

    DistRASCivector(const DistRASCivector<DataType>& o);
    DistRASCivector(std::shared_ptr<const DistRASCivector<DataType>> o) : DistRASCivector(*o) {}
    DistRASCivector(DistRASCivector<DataType>&& o);

    // Copy assignment
    DistRASCivector<DataType>& operator=(const DistRASCivector<DataType>& o);

    // Move assignment
    DistRASCivector<DataType>& operator=(DistRASCivector<DataType>&& o);

    using RASCivector_base<DistCIBlock<DataType>>::block;

    int get_bstring_buf(double* buf, const size_t a) const;

    void zero();

    void synchronize(const int root = 0) { /* do nothing */ }

    std::shared_ptr<DistRASCivector<DataType>> clone() const { return std::make_shared<DistRASCivector<DataType>>(det_); }
    std::shared_ptr<DistRASCivector<DataType>> copy() const  { return std::make_shared<DistRASCivector<DataType>>(*this); }
    std::shared_ptr<DistRASCivector<DataType>> transpose(std::shared_ptr<const RASDeterminants> det = nullptr) const;

//  TODO
//  std::shared_ptr<RASCivector<DataType>> civec() const { return std::make_shared<RASCivector<DataType>>(*this); }

    DataType dot_product(const DistRASCivector<DataType>& o) const;

    double norm() const { return std::sqrt(dot_product(*this)); }
    double variance() const { return dot_product(*this) / global_size_; }
    double rms() const { return std::sqrt(variance()); }

    void scale(const DataType a);
    void ax_plus_y(const DataType a, const DistRASCivector<DataType>& o);
    void ax_plus_y(const DataType a, std::shared_ptr<const DistRASCivector<DataType>> o) { ax_plus_y(a, *o); }

    // Spin functions are only implememted as specialized functions for double (see civec.cc)
    // returns < S^2 >
    DataType spin_expectation() const {
      std::shared_ptr<const DistRASCivector<DataType>> S2 = spin();
      return this->dot_product(*S2);
    }
    std::shared_ptr<DistRASCivector<DataType>> spin() const; // returns S^2 | civec >
    void spin_decontaminate(const double thresh = 1.0e-8);

    void project_out(std::shared_ptr<const DistRASCivector<DataType>> o) { ax_plus_y(-detail::conj(dot_product(*o)), *o); }

    double orthog(std::list<std::shared_ptr<const DistRASCivector<DataType>>> c) {
      for (auto& iter : c)
        project_out(iter);
      return normalize();
    }

    double orthog(std::shared_ptr<const DistRASCivector<DataType>> o) {
      return orthog(std::list<std::shared_ptr<const DistRASCivector<DataType>>>{o});
    }

    double normalize() {
      const double norm = this->norm();
      const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
      scale(DataType(scal));
      return norm;
    }

    void print(const double thr = 0.05) const;
};

template<> std::shared_ptr<DistRASCivector<double>> DistRASCivector<double>::spin() const; // returns S^2 | civec >
template<> void DistRASCivector<double>::spin_decontaminate(const double thresh);

using DistRASCivec = DistRASCivector<double>;
using DistRASDvec = Dvector_base<DistRASCivec>;

}

extern template class bagel::DistRASCivector<double>;

#endif
