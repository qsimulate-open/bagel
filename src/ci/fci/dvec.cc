//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dvec.cc
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

#include <stdexcept>
#include <src/ci/fci/dvec.h>

using namespace std;
using namespace bagel;

template <typename DataType>
Dvector<DataType>::Dvector(shared_ptr<const Determinants> det, const size_t ij)
  : btas::Tensor3<DataType>(det->lenb(), det->lena(), ij), det_(det), lena_(det->lena()), lenb_(det->lenb()), ij_(ij) {
  zero();
  DataType* tmp = data();
  for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_)
  dvec_.push_back(make_shared<Civector<DataType>>(det_, tmp));
}


template <typename DataType>
Dvector<DataType>::Dvector(const Dvector<DataType>& o) : btas::Tensor3<DataType>(o), det_(o.det_), lena_(o.lena_), lenb_(o.lenb_), ij_(o.ij_) {
  DataType* tmp = data();
  for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_)
    dvec_.push_back(make_shared<Civector<DataType>>(det_, tmp));
}


template <typename DataType>
Dvector<DataType>::Dvector(const Dvector_base<Civector<DataType>>& o)
 : btas::Tensor3<DataType>(o.det()->lenb(), o.det()->lena(), o.ij()), det_(o.det()), lena_(o.det()->lena()), lenb_(o.det()->lenb()), ij_(o.ij()) {
  DataType* tmp = data();
  for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_)
    dvec_.push_back(make_shared<Civector<DataType>>(det_, tmp));
  auto iter = dvec_.begin();
  for (auto& inp : o.dvec())
    **iter++ = *inp;
}


template <typename DataType>
Dvector<DataType>::Dvector(shared_ptr<const Civector<DataType>> e, const size_t ij)
  : btas::Tensor3<DataType>(e->lenb(), e->lena(), ij), det_(e->det()), lena_(e->lena()), lenb_(e->lenb()), ij_(ij) {
  DataType* tmp = data();
  for (int i = 0; i != ij; ++i, tmp += lenb_*lena_) {
    auto c = make_shared<Civector<DataType>>(det_, tmp);
    copy_n(e->data(), lenb_*lena_, c->data());
    dvec_.push_back(c);
  }
}


template <typename DataType>
vector<shared_ptr<const Civector<DataType>>> Dvector<DataType>::dvec(const vector<int>& conv) const {
  vector<shared_ptr<const Civector<DataType>>> out;
  auto c = conv.begin();
  for (auto& i : dvec_) {
    if (*c++ == 0) out.push_back(i);
    else out.push_back(nullptr);
  }
  return out;
}


template <typename DataType>
void Dvector<DataType>::set_det(shared_ptr<const Determinants> o) const {
  det_ = o;
  for_each(dvec_.begin(), dvec_.end(), [&o](CiPtr p){ p->set_det(o); });
}


template <typename DataType>
DataType Dvector<DataType>::dot_product(const Dvector<DataType>& other) const {
  return inner_product(dvec_.begin(), dvec_.end(), other.dvec_.begin(), DataType(0.0), plus<DataType>(), [](CiPtr p, CiPtr q){ return p->dot_product(q); });
}


template <typename DataType>
void Dvector<DataType>::ax_plus_y(const DataType a, const Dvector<DataType>& other) {
  transform(other.dvec_.begin(), other.dvec_.end(), dvec_.begin(), dvec_.begin(), [&a](CiPtr p, CiPtr q){ q->ax_plus_y(a, p); return q; });
}


template <typename DataType>
Dvector<DataType>& Dvector<DataType>::operator/=(const Dvector<DataType>& o) {
  assert(dvec().size() == o.dvec().size());
  transform(o.dvec_.begin(), o.dvec_.end(), dvec_.begin(), dvec_.begin(), [](CiPtr p, CiPtr q){ *q / *p; return q; });
  return *this;
}


template <typename DataType>
Dvector<DataType> Dvector<DataType>::operator/(const Dvector<DataType>& o) const {
  Dvector<DataType> out(*this);
  out /= o;
  return out;
}


template <typename DataType>
void Dvector<DataType>::scale(const DataType& a) {
  for_each(dvec_.begin(), dvec_.end(), [&a](CiPtr p){ p->scale(a); });
}


template <typename DataType>
void Dvector<DataType>::orthog(shared_ptr<const Dvector<DataType>> o) {
  if (o->ij() != ij())
    throw logic_error("Dvector<DataType>::orthog called inconsistently");
  transform(o->dvec_.begin(), o->dvec_.end(), dvec_.begin(), dvec_.begin(), [](CiPtr p, CiPtr q){ q->orthog(p); return q; });
}


template <typename DataType>
void Dvector<DataType>::project_out(shared_ptr<const Dvector<DataType>> o) {
  if (o->ij() != ij()) throw std::logic_error("Dvec::project_out called inconsistently");
  auto j = o->dvec().begin();
  for (auto& i : dvec())
    i->project_out(*j++);
}


template <typename DataType>
void Dvector<DataType>::project_out_all(shared_ptr<const Dvector<DataType>> o) {
  for (auto& i : dvec())
    for (auto& j : o->dvec())
      i->project_out(j);
}


template <typename DataType>
void Dvector<DataType>::synchronize() {
  for (auto& i : dvec_)
    i->synchronize();
}

template <typename DataType>
void Dvector<DataType>::rotate(std::shared_ptr<const Matrix> msrot) {
  // now coded in somewhat un-efficient manner
  Dvector<DataType> tmp(*this);
  const size_t detnumb = dvec().size();
  const size_t detsize = data(0)->size();
  for (size_t i = 0; i != detnumb; ++i) {
    for (size_t a = 0; a != detsize; ++a) {
      data(i)->data(a) = 0.0;
    }
  }
  for (size_t i = 0; i != detnumb; ++i) {
    for (size_t k = 0; k != detnumb; ++k) {
      for (size_t a = 0; a != detsize; ++a) {
        data(i)->data(a) += msrot->element(k, i) * tmp.data(k)->data(a);
      }
    }
  }
}


template <typename DataType>
void Dvector<DataType>::print(const double thresh) const {
  int j = 0;
  for (auto& iter : dvec_) {
    cout << endl << "     * ci vector, state " << setw(3) << j++;
    if (is_same<DataType, double>::value)
      cout << ", <S^2> = " << setw(6) << setprecision(4) << iter->spin_expectation();
    cout << endl;
    iter->print(thresh);
  }
}

template <typename DataType>
void Dvector<DataType>::print(const bool sort) const {
  int j = 0;
  for (auto& iter : dvec_) {
    std::cout << std::endl << "     * ci vector, state " << std::setw(3) << j++;
    if (typeid(DataType) == typeid(double)) std::cout << ", <S^2> = " << std::setw(6) << std::setprecision(4) << iter->spin_expectation();
    std::cout << std::endl;
    iter->print(0.0, false);
  }
}


template <typename DataType>
shared_ptr<Dvector<DataType>> Dvector<DataType>::extract_state(const vector<int> input) const {
  auto out = make_shared<Dvector<DataType>>(det(), input.size());
  for (int i = 0; i != input.size(); ++i)
    copy_n(data(input[i])->data(), lenb_*lena_, out->data(i)->data());
    //out->data(i) = data(input[i])->copy();
  return out;
}


template class bagel::Dvector<double>;
template class bagel::Dvector<complex<double>>;

