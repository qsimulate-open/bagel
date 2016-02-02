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
#ifdef HAVE_GA
#include <ga.h>
#endif
#include <src/ci/fci/civec.h>

namespace bagel {

// task class
template<typename DataType>
class GA_Task {
#ifndef HAVE_GA
  public:
    using ga_nbhdl_t = long; // jsut to compile
#endif
  protected:
    ga_nbhdl_t tag;
    std::unique_ptr<DataType[]> buf;
  public:
    GA_Task(ga_nbhdl_t t) : tag(t) { }
    GA_Task(ga_nbhdl_t t, std::unique_ptr<DataType[]>&& o) : tag(t), buf(std::move(o)) { }
    void wait();
};

extern template class GA_Task<double>;
extern template class GA_Task<std::complex<double>>;


template<typename DataType>
class DistCivector {
  public: using DetType = Determinants;
  public: using LocalizedType = std::false_type;

  protected:
    mutable std::shared_ptr<const Determinants> det_;

    // global dimension
    size_t lena_;
    size_t lenb_;

    // global array tag
    int ga_;

    // table for alpha string distribution
    StaticDist dist_;

    // local alpha strings
    size_t astart_;
    size_t aend_;

    // distribution information
    std::vector<int64_t> blocks_;

  public:
    DistCivector(std::shared_ptr<const Determinants> det);

    DistCivector(const DistCivector<DataType>& o) : DistCivector(o.det()) { *this = o; }
    DistCivector(std::shared_ptr<const DistCivector<DataType>> o) : DistCivector(*o) {}

    DistCivector<DataType>& operator=(const DistCivector<DataType>& o);

    std::shared_ptr<DistCivector<DataType>> clone() const { return std::make_shared<DistCivector<DataType>>(det_); }
    std::shared_ptr<DistCivector<DataType>> copy() const { return std::make_shared<DistCivector<DataType>>(*this); }

    size_t size() const { return lenb_*(aend_-astart_); }
    size_t global_size() const { return lena_*lenb_; }
    size_t lena() const { return lena_; }
    size_t lenb() const { return lenb_; }

    size_t astart() const { return astart_; }
    size_t aend() const { return aend_; }
    size_t asize() const { return aend_ - astart_; }

    void zero();
    void synchronize(const int root = 0) { /* do nothing */ }

    std::shared_ptr<Civector<DataType>> civec() const;
    std::shared_ptr<const Determinants> det() const { return det_; }

    void set_det(std::shared_ptr<const Determinants> d) { det_ = d; }

    bool is_local(const size_t a) const;
    std::unique_ptr<DataType[]> local() const;
    void set_local(const size_t la, const size_t lb, const DataType a);

    std::shared_ptr<GA_Task<DataType>> accumulate_bstring_buf(std::unique_ptr<DataType[]>&& buf, const size_t a);
    std::shared_ptr<GA_Task<DataType>> get_bstring_buf(DataType* buf, const size_t a) const;

    void local_accumulate(const DataType a, const std::unique_ptr<DataType[]>& buf);

    // utility functions
    DataType dot_product(const DistCivector<DataType>& o) const;
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    double variance() const { return detail::real(dot_product(*this)) / (lena_*lenb_); }
    double rms() const { return std::sqrt(variance()); }

    void scale(const DataType a);
    void ax_plus_y(const DataType a, const DistCivector<DataType>& o);
    void ax_plus_y(const DataType a, std::shared_ptr<const DistCivector<DataType>> o) { ax_plus_y(a, *o); }
    void project_out(std::shared_ptr<const DistCivector<DataType>> o) { ax_plus_y(-detail::conj(dot_product(*o)), *o); }

#if 0
    DataType spin_expectation() const {
      std::shared_ptr<DistCivector<DataType>> S2 = spin();
      return dot_product(*S2);
    }
    std::shared_ptr<DistCivector<DataType>> spin() const;
    void spin_decontaminate(const double thresh = 1.0e-4);
#else
    DataType spin_expectation() const { return 0.0; }
    std::shared_ptr<DistCivector<DataType>> spin() const { return nullptr; }
    void spin_decontaminate(const double thresh = 1.0e-4) { }
#endif

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

#if 0
template <> std::shared_ptr<DistCivector<double>> DistCivector<double>::spin() const;
template <> void DistCivector<double>::spin_decontaminate(const double);
#endif

template <> double DistCivector<double>::dot_product(const DistCivector<double>&) const;
template <> std::complex<double> DistCivector<std::complex<double>>::dot_product(const DistCivector<std::complex<double>>&) const;

extern template class DistCivector<double>;
extern template class DistCivector<std::complex<double>>;

using DistCivec = DistCivector<double>;
using ZDistCivec = DistCivector<std::complex<double>>;

using DistDvec = Dvector_base<DistCivec>;

}

#endif
