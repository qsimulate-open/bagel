//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: civec.h
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


#ifndef BAGEL_FCI_CIVEC_H
#define BAGEL_FCI_CIVEC_H

#include <list>
#include <src/util/math/algo.h>
#include <src/util/f77.h>
#include <src/util/parallel/staticdist.h>
#include <src/ci/fci/determinants.h>
#include <src/ci/fci/dvector_base.h>

namespace bagel {

template<typename DataType>
class Civector {
  public: using DetType = Determinants; // used to automatically determine type for Determinants object in templates
  public: using LocalizedType = std::true_type;

  protected:
    // The determinant space in which this Civec object is defined
    mutable std::shared_ptr<const Determinants> det_;

    size_t lena_;
    size_t lenb_;

    // !!CAUTION!!
    // cc is formated so that B runs first.
    // Also, cc_ can be null if this is constructed by Dvec.
    std::unique_ptr<DataType[]> cc_;

    DataType* cc_ptr_;

    DataType& cc(const size_t& i) { return *(cc_ptr_+i); }
    const DataType& cc(const size_t& i) const { return *(cc_ptr_+i); }
    DataType* cc() { return cc_ptr_; }
    const DataType* cc() const { return cc_ptr_; }

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      if (!cc_.get())
        throw std::logic_error("illegal call of Civector<T>::save");
      ar << det_ << lena_ << lenb_ << make_array(cc(), size());
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> det_ >> lena_ >> lenb_;
      cc_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      cc_ptr_ = cc_.get();
      ar >> make_array(cc(), size());
    }

  public:
    Civector() { }

    Civector(std::shared_ptr<const Determinants> det) : det_(det), lena_(det->lena()), lenb_(det->lenb()) {
      cc_ = std::unique_ptr<DataType[]>(new DataType[lena_*lenb_]);
      cc_ptr_ = cc_.get();
      std::fill_n(cc(), lena_*lenb_, 0.0);
    }

    // constructor that is called by Dvec.
    Civector(std::shared_ptr<const Determinants> det, DataType* din_) : det_(det), lena_(det->lena()), lenb_(det->lenb()) {
      cc_ptr_ = din_;
    }

    // copy constructor
    Civector(const Civector<DataType>& o) : det_(o.det_), lena_(o.lena_), lenb_(o.lenb_) {
      cc_ = std::unique_ptr<DataType[]>(new DataType[lena_*lenb_]);
      cc_ptr_ = cc_.get();
      std::copy_n(o.cc(), lena_*lenb_, cc());
    }

    // this is not a copy constructor.
    Civector(std::shared_ptr<Civector<DataType>> o, std::shared_ptr<const Determinants> det) : det_(det), lena_(o->lena_), lenb_(o->lenb_) {
      assert(lena_ == det->lena() && lenb_ == det->lenb());
      cc_ = std::move(o->cc_);
      cc_ptr_ = cc_.get();
    }

    std::shared_ptr<Civector<DataType>> clone() const { return std::make_shared<Civector<DataType>>(det_); }
    std::shared_ptr<Civector<DataType>> copy() const { return std::make_shared<Civector<DataType>>(*this); }

    DataType* data() { return cc(); }
    DataType& element(size_t i, size_t j) { return cc(i+j*lenb_); }
    DataType* element_ptr(size_t i, size_t j) { return cc()+i+j*lenb_; }

    DataType& data(const size_t& i) { return cc(i); }
    const DataType& data(const size_t& i) const { return cc(i); }

    const DataType* data() const { return cc(); }
    const DataType* element_ptr(size_t i, size_t j) const { return cc()+i+j*lenb_; }

    std::shared_ptr<const Determinants> det() const { return det_; }
    void set_det(std::shared_ptr<const Determinants> o) const { det_ = o; }

    void zero() { std::fill(cc(), cc()+lena_*lenb_, DataType(0.0)); }

    size_t size() const { return lena_*lenb_; }

    std::shared_ptr<Civector<DataType>> transpose(std::shared_ptr<const Determinants> det = nullptr) const {
      if (det == nullptr) det = det_->transpose();
      auto ct = std::make_shared<Civector<DataType>>(det);
      blas::transpose(cc(), lenb_, lena_, ct->data());

      if (det_->nelea()*det_->neleb() & 1)
        ct->scale(-1.0);
      return ct;
    }

    size_t lena() const { return lena_; }
    size_t lenb() const { return lenb_; }

    // some functions for convenience
    void ax_plus_y(DataType a, const Civector<DataType>& other) {
      assert((lena_ == other.lena_) && (lenb_ == other.lenb_));
      std::transform(other.data(), other.data()+size(), cc(), cc(), [&a](DataType p, DataType q){ return a*p+q; });
    }
    void ax_plus_y(DataType a, std::shared_ptr<const Civector> other) { ax_plus_y(a, *other); }

    DataType dot_product(const Civector<DataType>& other) const {
      assert((lena_ == other.lena_) && (lenb_ == other.lenb_));
      return blas::dot_product(cc(), size(), other.data());
    }
    DataType dot_product(std::shared_ptr<const Civector> other) const { return dot_product(*other); }

    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    double variance() const { return detail::real(dot_product(*this)) / size(); }
    double rms() const { return std::sqrt(variance()); }

    void scale(const DataType a) {
      std::for_each(cc(), cc()+size(), [&a](DataType& p){ p *= a; });
    }

    // Spin functions are only implememted as specialized functions for double (see civec.cc)
    DataType spin_expectation() const { // returns < S^2 >
      std::shared_ptr<Civector<DataType>> S2 = spin();
      return dot_product(*S2);
    }
    std::shared_ptr<Civector<DataType>> spin() const {
      auto out = std::make_shared<Civector<DataType>>(det_);

      // First the easy part, S_z^2 + S_z
      const double sz = 0.5*static_cast<double>(det_->nspin());
      *out = *this;
      *out *= sz*sz + sz + det_->neleb();

      const int norb = det_->norb();
      const int lena = det_->lena();
      const int lenb = det_->lenb();

      auto intermediate = std::make_shared<Civector<DataType>>(det_);

      for (int i = 0; i < norb; ++i) {
        for (int j = 0; j < norb; ++j) {
          intermediate->zero();
          for (auto& iter : det_->phia(i,j)) {
            const DataType* source = this->element_ptr(0, iter.source);
            DataType* target = intermediate->element_ptr(0, iter.target);
            double sign = static_cast<double>(iter.sign);

            std::transform(source, source+lenb, target, target, [&sign](DataType p, DataType q){ return q+sign*p; });
          }
          for (int ia = 0; ia < lena; ++ia) {
            DataType* target_base = out->element_ptr(0, ia);
            const DataType* source_base = intermediate->element_ptr(0, ia);
            for (auto& iter : det_->phib(j,i)) {
              target_base[iter.target] -= static_cast<double>(iter.sign) * source_base[iter.source];
            }
          }
        }
      }
      return out;
    }

    // S_- = \sum_i i_beta^\dagger i_alpha
    std::shared_ptr<Civector<DataType>> spin_lower(std::shared_ptr<const Determinants> target_det = nullptr) const {
      if (target_det == nullptr)
        target_det = std::make_shared<Determinants>(det_->norb(), det_->nelea()-1, det_->neleb()+1, det_->compress(), true);
      assert( (target_det->nelea() == det_->nelea()-1) && (target_det->neleb() == det_->neleb()+1) );
      auto out = std::make_shared<Civector<DataType>>(target_det);
      std::shared_ptr<const Determinants> source_det = det_;
      const int norb = source_det->norb();
      const int source_lena = source_det->lena();
      const int source_lenb = source_det->lenb();

      DataType* source_data = cc_ptr_;
      // This is a safe but probably slow implementation
      for (int aiter = 0; aiter < source_lena; ++aiter) {
        const std::bitset<nbit__> alphastring = source_det->string_bits_a(aiter);
        for (int biter = 0; biter < source_lenb; ++biter, ++source_data) {
          const std::bitset<nbit__> betastring = source_det->string_bits_b(biter);
          for (int i = 0; i < norb; ++i) {
            std::bitset<nbit__> abit = alphastring;
            std::bitset<nbit__> bbit = betastring;
            if (abit[i]) {
              abit.reset(i);
              if (!bbit[i]) {
                bbit.set(i);

                const int atarget = target_det->lexical<0>(abit);
                const int btarget = target_det->lexical<1>(bbit);
                // Now the computation begins

                const int aphase = source_det->sign<0>(alphastring, i);
                const int bphase = -1*source_det->sign<1>(betastring, i);

                out->element(btarget, atarget) += static_cast<double>(aphase*bphase) * (*source_data);
              }
            }
          }
        }
      }
      return out;
    }

    // S_+ = \sum_i i_alpha^\dagger i_beta
    std::shared_ptr<Civector<DataType>> spin_raise(std::shared_ptr<const Determinants> target_det = nullptr) const {
      if (target_det == nullptr)
        target_det = std::make_shared<Determinants>(det_->norb(), det_->nelea()+1, det_->neleb()-1, det_->compress(), true);
      assert( (target_det->nelea() == det_->nelea()+1) && (target_det->neleb() == det_->neleb()-1) );
      auto out = std::make_shared<Civector<DataType>>(target_det);

      std::shared_ptr<const Determinants> source_det = det_;
      const int norb = source_det->norb();
      const int source_lena = source_det->lena();
      const int source_lenb = source_det->lenb();

      DataType* source_data = cc_ptr_;
      // This is a safe but probably slow implementation
      for (int aiter = 0; aiter < source_lena; ++aiter) {
        const std::bitset<nbit__> alphastring = source_det->string_bits_a(aiter);
        for (int biter = 0; biter < source_lenb; ++biter, ++source_data) {
          const std::bitset<nbit__> betastring = source_det->string_bits_b(biter);
          for (int i = 0; i < norb; ++i) {
            std::bitset<nbit__> abit = alphastring;
            std::bitset<nbit__> bbit = betastring;
            if (bbit[i]) {
              bbit.reset(i);
              if (!abit[i]) {
                abit.set(i);

                const int atarget = target_det->lexical<0>(abit);
                const int btarget = target_det->lexical<1>(bbit);

                const int aphase = source_det->sign<0>(alphastring, i);
                const int bphase = source_det->sign<1>(betastring, i);

                out->element(btarget, atarget) += static_cast<double>(aphase*bphase) * (*source_data);
              }
            }
          }
        }
      }
      return out;
    }

    void spin_decontaminate(const double thresh = 1.0e-12);

    std::shared_ptr<Civector<DataType>> apply(const int orbital, const bool action, const bool spin) const {
      // action: true -> create; false -> annihilate
      // spin: true -> alpha; false -> beta

      std::shared_ptr<const Determinants> source_det = this->det();

      const int source_lenb = source_det->lenb();

      std::shared_ptr<Civector<DataType>> out;

      if (spin) {
        std::shared_ptr<const Determinants> target_det = ( action ? source_det->addalpha() : source_det->remalpha() );
        out = std::make_shared<Civector<DataType>>(target_det);

        const int target_lenb = target_det->lenb();

        DataType* target_base = out->data();
        const DataType* source_base = this->data();
        for (auto& iter : ( action ? source_det->phiupa(orbital) : source_det->phidowna(orbital) )) {
          const DataType sign = static_cast<DataType>(iter.sign);
          DataType* target = target_base + target_lenb * iter.target;
          const DataType* source = source_base + source_lenb * iter.source;
          std::transform(source, source + target_lenb, target, target, [&sign] (DataType p, DataType q) { return sign * p + q; });
        }
      }
      else {
        std::shared_ptr<const Determinants> target_det = ( action ? source_det->addbeta() : source_det->rembeta() );

        const int target_lena = target_det->lena();

        out = std::make_shared<Civector<DataType>>(target_det);

        for (int i = 0; i < target_lena; ++i) {
          DataType* target_base = out->element_ptr(0,i);
          const DataType* source_base = this->element_ptr(0,i);
          for (auto& iter : ( action ? source_det->phiupb(orbital) : source_det->phidownb(orbital) )) {
            const DataType sign = static_cast<DataType>(iter.sign);
            target_base[iter.target] += sign * source_base[iter.source];
          }
        }
      }

      return out;
    }

    Civector<DataType>& operator*=(const double& a) { scale(a); return *this; }
    Civector<DataType>& operator+=(const double& a) { std::for_each(cc(), cc()+size(), [&a](DataType& p){ p += a; }); return *this; }
    Civector<DataType>& operator-=(const double& a) { std::for_each(cc(), cc()+size(), [&a](DataType& p){ p -= a; }); return *this; }

    Civector<DataType>& operator=(const Civector<DataType>& o) { assert(det()->lena() == o.det()->lena() && det()->lenb() == o.det()->lenb()); std::copy_n(o.cc(), size(), cc()); return *this; }
    Civector<DataType>& operator+=(const Civector<DataType>& o) { ax_plus_y( 1.0, o); return *this; }
    Civector<DataType>& operator-=(const Civector<DataType>& o) { ax_plus_y(-1.0, o); return *this; }
    Civector<DataType>& operator/=(const Civector<DataType>& o) {
      for (size_t i = 0; i != size(); ++i)
        data(i) /= o.data(i);
      return *this;
    }
    Civector<DataType> operator/(const Civector<DataType>& o) const {
      Civector<DataType> out(*this);
      out /= o;
      return out;
    }

    // assumes that Civec's in c are already orthogonal with each other.
    // returns scaling factor (see implementation)

    double orthog(std::list<std::shared_ptr<const Civector<DataType>>> c) {
      for (auto& iter : c)
        project_out(iter);
      return normalize();
    }

    double orthog(std::shared_ptr<const Civector<DataType>> o) {
      return orthog(std::list<std::shared_ptr<const Civector<DataType>>>{o});
    }

    double normalize() {
      const double norm = this->norm();
      const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
      scale(scal);
      return 1.0/scal;
    }

    void project_out(std::shared_ptr<const Civector<DataType>> o) { ax_plus_y(-detail::conj(dot_product(*o)), *o); }

    void print(const double thr) const {
      const DataType* i = cc();
      // multimap sorts elements so that they will be shown in the descending order in magnitude
      std::multimap<double, std::tuple<DataType, std::bitset<nbit__>, std::bitset<nbit__>>> tmp;
      for (auto& ia : det_->string_bits_a()) {
        for (auto& ib : det_->string_bits_b()) {
          if (std::abs(*i) > thr)
            tmp.emplace(-std::abs(*i), std::make_tuple(*i, ia, ib));
          ++i;
        }
      }
      for (auto& iter : tmp)
        std::cout << "       " << print_bit(std::get<1>(iter.second), std::get<2>(iter.second), det()->norb())
                  << "  " << std::setprecision(10) << std::setw(15) << std::get<0>(iter.second) << std::endl;
    }

    void synchronize(const int root = 0) {
#ifdef HAVE_MPI_H
      mpi__->broadcast(cc_ptr_, size(), root);
#endif
    }
};

template<> void Civector<double>::spin_decontaminate(const double thresh);

using Civec = Civector<double>;
using ZCivec = Civector<std::complex<double>>;

template <>
void Dvector_base<Civec>::apply_and_fill(std::shared_ptr<const Dvector_base<Civec>> source_dvec, const int orbital, const bool action, const bool spin);

using CASDvec = Dvector_base<Civec>;

}

extern template class bagel::Civector<double>;
extern template class bagel::Civector<std::complex<double>>;

#endif
