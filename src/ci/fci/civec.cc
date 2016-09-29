//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: civec.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/ci/fci/civec.h>
#include <src/util/math/algo.h>

using namespace std;
using namespace bagel;

template<typename DataType>
Civector<DataType>::Civector(shared_ptr<const Determinants> det) : det_(det), lena_(det->lena()), lenb_(det->lenb()) {
  cc_ = unique_ptr<DataType[]>(new DataType[lena_*lenb_]);
  cc_ptr_ = cc_.get();
  fill_n(cc(), lena_*lenb_, 0.0);
}


// constructor that is called by Dvec.
template<typename DataType>
Civector<DataType>::Civector(shared_ptr<const Determinants> det, DataType* din_) : det_(det), lena_(det->lena()), lenb_(det->lenb()) {
  cc_ptr_ = din_;
}


// copy constructor
template<typename DataType>
Civector<DataType>::Civector(const Civector<DataType>& o) : det_(o.det_), lena_(o.lena_), lenb_(o.lenb_) {
  cc_ = unique_ptr<DataType[]>(new DataType[lena_*lenb_]);
  cc_ptr_ = cc_.get();
  copy_n(o.cc(), lena_*lenb_, cc());
}


// this is not a copy constructor.
template<typename DataType>
Civector<DataType>::Civector(shared_ptr<Civector<DataType>> o, shared_ptr<const Determinants> det) : det_(det), lena_(o->lena_), lenb_(o->lenb_) {
  assert(lena_ == det->lena() && lenb_ == det->lenb());
  cc_ = move(o->cc_);
  cc_ptr_ = cc_.get();
}


template<typename DataType>
shared_ptr<Civector<DataType>> Civector<DataType>::transpose(shared_ptr<const Determinants> det) const {
  if (det == nullptr) det = det_->transpose();
  auto ct = make_shared<Civector<DataType>>(det);
  blas::transpose(cc(), lenb_, lena_, ct->data());

  if (det_->nelea()*det_->neleb() & 1)
    ct->scale(-1.0);
  return ct;
}


template<typename DataType>
shared_ptr<Civector<DataType>> Civector<DataType>::spin() const {
  auto out = make_shared<Civector<DataType>>(det_);

  // First the easy part, S_z^2 + S_z
  const double sz = 0.5*static_cast<double>(det_->nspin());
  *out = *this;
  *out *= sz*sz + sz + det_->neleb();

  const int norb = det_->norb();
  const int lena = det_->lena();
  const int lenb = det_->lenb();

  auto intermediate = make_shared<Civector<DataType>>(det_);

  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      intermediate->zero();
      for (auto& iter : det_->phia(i,j)) {
        const DataType* source = this->element_ptr(0, iter.source);
        DataType* target = intermediate->element_ptr(0, iter.target);
        blas::ax_plus_y_n(static_cast<double>(iter.sign), source, lenb, target);
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


template<typename DataType>
shared_ptr<Civector<DataType>> Civector<DataType>::spin_lower(shared_ptr<const Determinants> target_det) const {
  if (target_det == nullptr)
    target_det = make_shared<Determinants>(det_->norb(), det_->nelea()-1, det_->neleb()+1, det_->compress(), true);
  assert( (target_det->nelea() == det_->nelea()-1) && (target_det->neleb() == det_->neleb()+1) );
  auto out = make_shared<Civector<DataType>>(target_det);
  shared_ptr<const Determinants> source_det = det_;
  const int norb = source_det->norb();
  const int source_lena = source_det->lena();
  const int source_lenb = source_det->lenb();

  DataType* source_data = cc_ptr_;
  // This is a safe but probably slow implementation
  for (int aiter = 0; aiter < source_lena; ++aiter) {
    const bitset<nbit__> alphastring = source_det->string_bits_a(aiter);
    for (int biter = 0; biter < source_lenb; ++biter, ++source_data) {
      const bitset<nbit__> betastring = source_det->string_bits_b(biter);
      for (int i = 0; i < norb; ++i) {
        bitset<nbit__> abit = alphastring;
        bitset<nbit__> bbit = betastring;
        if (abit[i] && !bbit[i]) {
          abit.reset(i);
          bbit.set(i);
          const int atarget = target_det->lexical<0>(abit);
          const int btarget = target_det->lexical<1>(bbit);
          const double abphase = -1*source_det->sign<0>(alphastring, i)*source_det->sign<1>(betastring, i);
          out->element(btarget, atarget) += abphase * *source_data;
        }
      }
    }
  }
  return out;
}


template<typename DataType>
shared_ptr<Civector<DataType>> Civector<DataType>::spin_raise(shared_ptr<const Determinants> target_det) const {
  if (target_det == nullptr)
    target_det = make_shared<Determinants>(det_->norb(), det_->nelea()+1, det_->neleb()-1, det_->compress(), true);
  assert(target_det->nelea() == det_->nelea()+1 && target_det->neleb() == det_->neleb()-1);
  auto out = make_shared<Civector<DataType>>(target_det);

  shared_ptr<const Determinants> source_det = det_;
  const int norb = source_det->norb();
  const int source_lena = source_det->lena();
  const int source_lenb = source_det->lenb();

  DataType* source_data = cc_ptr_;
  // This is a safe but probably slow implementation
  for (int aiter = 0; aiter < source_lena; ++aiter) {
    const bitset<nbit__> alphastring = source_det->string_bits_a(aiter);
    for (int biter = 0; biter < source_lenb; ++biter, ++source_data) {
      const bitset<nbit__> betastring = source_det->string_bits_b(biter);
      for (int i = 0; i < norb; ++i) {
        bitset<nbit__> abit = alphastring;
        bitset<nbit__> bbit = betastring;
        if (bbit[i] && !abit[i]) {
          abit.set(i);
          bbit.reset(i);
          const int atarget = target_det->lexical<0>(abit);
          const int btarget = target_det->lexical<1>(bbit);
          const double abphase = source_det->sign<0>(alphastring, i) * source_det->sign<1>(betastring, i);
          out->element(btarget, atarget) += abphase * *source_data;
        }
      }
    }
  }
  return out;
}


template<typename DataType>
shared_ptr<Civector<DataType>> Civector<DataType>::apply(const int orbital, const bool action, const bool spin) const {
  // action: true -> create; false -> annihilate
  // spin: true -> alpha; false -> beta
  shared_ptr<const Determinants> source_det = this->det();
  const int source_lenb = source_det->lenb();

  shared_ptr<Civector<DataType>> out;
  if (spin) {
    shared_ptr<const Determinants> target_det = action ? source_det->addalpha() : source_det->remalpha();
    out = make_shared<Civector<DataType>>(target_det);
    const int target_lenb = target_det->lenb();

    DataType* target_base = out->data();
    const DataType* source_base = this->data();
    for (auto& iter : (action ? source_det->phiupa(orbital) : source_det->phidowna(orbital))) {
      const DataType sign = static_cast<DataType>(iter.sign);
      DataType* target = target_base + target_lenb * iter.target;
      const DataType* source = source_base + source_lenb * iter.source;
      blas::ax_plus_y_n(sign, source, target_lenb, target);
    }
  } else {
    shared_ptr<const Determinants> target_det = action ? source_det->addbeta() : source_det->rembeta();
    out = make_shared<Civector<DataType>>(target_det);
    const int target_lena = target_det->lena();

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


template<>
void Civector<double>::spin_decontaminate(const double thresh) {
  const int nspin = det_->nspin();
  const int max_spin = det_->nelea() + det_->neleb();
  const double expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<Civec> S2 = spin();

  int k = nspin + 2;
  while(fabs(dot_product(*S2) - expectation) > thresh) {
    if (k > max_spin) throw runtime_error("Spin decontamination failed.");
    const double factor = -4.0/(static_cast<double>(k*(k+2)));
    ax_plus_y(factor, *S2);
    const double norm = this->norm();
    const double rescale = (norm*norm > 1.0e-60) ? 1.0/norm : 0.0;
    scale(rescale);

    S2 = spin();
    k += 2;
  }
}


template<>
void Civector<complex<double>>::spin_decontaminate(const double thresh) {
  assert(false);
}


template<typename DataType>
void Civector<DataType>::print(const double thr, const bool sort) const {
  const DataType* i = cc();
  // multimap sorts elements so that they will be shown in the descending order in magnitude
  if (sort) {
    multimap<double, tuple<DataType, bitset<nbit__>, bitset<nbit__>>> tmp;
    for (auto& ia : det_->string_bits_a()) {
      for (auto& ib : det_->string_bits_b()) {
        if (abs(*i) > thr)
          tmp.emplace(-abs(*i), make_tuple(*i, ia, ib));
        ++i;
      }
    }
    for (auto& iter : tmp)
      cout << "       " << print_bit(get<1>(iter.second), get<2>(iter.second), det()->norb())
           << "  " << setprecision(10) << setw(15) << get<0>(iter.second) << endl;
  } else {
    for (auto& ia : det_->string_bits_a()) {
      for (auto& ib : det_->string_bits_b()) {
        cout << "       " << print_bit(ia, ib, det()->norb())
             << "  " << setprecision(10) << setw(15) << *i << endl;
        ++i;
      }
    }
  }
}


template<>
void CASDvec::apply_and_fill(shared_ptr<const CASDvec> source_dvec, const int orbital, const bool action, const bool spin) {
  shared_ptr<const DetType> source_det = source_dvec->det();
  shared_ptr<const DetType> target_det = this->det();

  this->zero();

  const int source_lenb = source_det->lenb();
  const int target_lenb = target_det->lenb();
  const int target_lena = target_det->lena();

  if (spin) {
    for (size_t ivec = 0; ivec != this->ij(); ++ivec) {
      double* target_base = this->data(ivec)->data();
      const double* source_base = source_dvec->data(ivec)->data();
      for (auto& iter : (action ? source_det->phiupa(orbital) : source_det->phidowna(orbital))) {
        const double sign = static_cast<double>(iter.sign);
        double* target = target_base + target_lenb * iter.target;
        const double* source = source_base + source_lenb * iter.source;
        blas::ax_plus_y_n(sign, source, target_lenb, target);
      }
    }
  } else {
    for (size_t ivec = 0; ivec != this->ij(); ++ivec) {
      for (int i = 0; i < target_lena; ++i) {
        double* target_base = this->data(ivec)->element_ptr(0,i);
        const double* source_base = source_dvec->data(ivec)->element_ptr(0,i);
        for (auto& iter : (action ? source_det->phiupb(orbital) : source_det->phidownb(orbital))) {
          const double sign = static_cast<double>(iter.sign);
          target_base[iter.target] += sign * source_base[iter.source];
        }
      }
    }
  }
}

template class bagel::Civector<double>;
template class bagel::Civector<complex<double>>;

