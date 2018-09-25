//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rotfile.cc
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

#include <src/multi/casscf/rotfile.h>

using namespace std;
using namespace bagel;


template<typename DataType>
RotationMatrix<DataType>::RotationMatrix(const int iclos, const int iact, const int ivirt)
  : nclosed_(iclos), nact_(iact), nvirt_(ivirt), size_(iclos*iact+iclos*ivirt+iact*ivirt), data_(new DataType[size_]) {
  zero();
}


template<typename DataType>
RotationMatrix<DataType>::RotationMatrix(const RotationMatrix& o)
  : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), size_(o.size_), data_(new DataType[o.size_]) {
  *this = o;
}


template<typename DataType>
RotationMatrix<DataType>::RotationMatrix(shared_ptr<const RotationMatrix> o)
  : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), size_(o->size_), data_(new DataType[o->size_]) {
  *this = *o;
}


template<typename DataType>
RotationMatrix<DataType>::RotationMatrix(shared_ptr<const Matrix_base<DataType>> o, const int iclos, const int iact, const int ivirt)
  : nclosed_(iclos), nact_(iact), nvirt_(ivirt), size_(iclos*iact+iclos*ivirt+iact*ivirt), data_(new DataType[size_]) {
  const int nocc = nclosed_ + nact_;
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      ele_va(j, i) = o->element(j+nocc, i+nclosed_);
    }
    for (int j = 0; j != nclosed_; ++j) {
      ele_ca(j, i) = o->element(i+nclosed_, j);
    }
  }
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      ele_vc(j, i) = o->element(j+nocc, i);
    }
  }
}


template<typename DataType>
shared_ptr<RotationMatrix<DataType>> RotationMatrix<DataType>::clone() const {
  return make_shared<RotationMatrix<DataType>>(nclosed_, nact_, nvirt_);
}


template<typename DataType>
shared_ptr<RotationMatrix<DataType>> RotationMatrix<DataType>::copy() const {
  return make_shared<RotationMatrix<DataType>>(*this);
}


// overloaded operators
template<typename DataType>
RotationMatrix<DataType> RotationMatrix<DataType>::operator+(const RotationMatrix<DataType>& o) const {
  RotationMatrix<DataType> out(*this);
  out.ax_plus_y(1.0, o);
  return out;
}


template<typename DataType>
RotationMatrix<DataType> RotationMatrix<DataType>::operator-(const RotationMatrix<DataType>& o) const {
  RotationMatrix<DataType> out(*this);
  out.ax_plus_y(-1.0, o);
  return out;
}


template<typename DataType>
RotationMatrix<DataType>& RotationMatrix<DataType>::operator+=(const RotationMatrix<DataType>& o) {
  ax_plus_y(1.0, o);
  return *this;
}


template<typename DataType>
RotationMatrix<DataType>& RotationMatrix<DataType>::operator-=(const RotationMatrix<DataType>& o) {
  ax_plus_y(-1.0, o);
  return *this;
}


template<typename DataType>
RotationMatrix<DataType>& RotationMatrix<DataType>::operator*=(const DataType a) {
  scale(a);
  return *this;
}


template<typename DataType>
RotationMatrix<DataType>& RotationMatrix<DataType>::operator/=(const RotationMatrix<DataType>& o) {
  for (int i = 0; i != size(); ++i)
    data(i)/= o.data(i);
  return *this;
}


template<typename DataType>
RotationMatrix<DataType>& RotationMatrix<DataType>::operator*=(const RotationMatrix<DataType>& o) {
  for (int i = 0; i != size(); ++i)
    data(i)*= o.data(i);
  return *this;
}


template<typename DataType>
RotationMatrix<DataType> RotationMatrix<DataType>::operator/(const RotationMatrix<DataType>& o) const {
  RotationMatrix<DataType> out(*this);
  return out /= o;
}


template<typename DataType>
RotationMatrix<DataType> RotationMatrix<DataType>::operator*(const RotationMatrix<DataType>& o) const {
  RotationMatrix<DataType> out(*this);
  return out *= o;
}


template<typename DataType>
RotationMatrix<DataType>& RotationMatrix<DataType>::operator=(const RotationMatrix<DataType>& o) {
  copy_n(o.data(), size(), data());
  return *this;
}


// orthogonalize to the list of RotationMatrix's
template<typename DataType>
double RotationMatrix<DataType>::orthog(list<shared_ptr<const RotationMatrix<DataType>>> c) {
  for (auto iter = c.begin(); iter != c.end(); ++iter)
    this->ax_plus_y(- detail::conj(this->dot_product(**iter)), **iter);
  return normalize();
}


template<typename DataType>
double RotationMatrix<DataType>::normalize() {
  const double scal = 1.0/this->norm();
  scale(scal);
  return 1.0/scal;
}


template<typename DataType>
void RotationMatrix<DataType>::synchronize() {
#ifdef HAVE_MPI_H
  mpi__->broadcast(data(), size(), 0);
#endif
}


template<typename DataType>
void RotationMatrix<DataType>::ax_plus_y_ca(const DataType a, const ViewType mat) {
  assert(mat.ndim() == nclosed_ && mat.mdim() == nact_);
  blas::ax_plus_y_n(a, mat.data(), nclosed_*nact_, ptr_ca());
}


template<typename DataType>
void RotationMatrix<DataType>::ax_plus_y_va(const DataType a, const ViewType mat) {
  assert(mat.ndim() == nvirt_ && mat.mdim() == nact_);
  blas::ax_plus_y_n(a, mat.data(), nvirt_*nact_,  ptr_va());
}


template<typename DataType>
void RotationMatrix<DataType>::ax_plus_y_vc(const DataType a, const ViewType mat) {
  assert(mat.ndim() == nvirt_ && mat.mdim() == nclosed_);
  blas::ax_plus_y_n(a, mat.data(), nvirt_*nclosed_, ptr_vc());
}


template<typename DataType>
shared_ptr<typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type>
  RotationMatrix<DataType>::ca_mat() const {
  auto out = make_shared<MatType>(nclosed_, nact_);
  copy_n(ptr_ca(), nclosed_*nact_, out->data());
  return out;
}


template<typename DataType>
shared_ptr<typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type>
  RotationMatrix<DataType>::va_mat() const {
  auto out = make_shared<MatType>(nvirt_, nact_);
  copy_n(ptr_va(), nvirt_*nact_, out->data());
  return out;
}


template<typename DataType>
shared_ptr<typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type>
  RotationMatrix<DataType>::vc_mat() const {
  auto out = make_shared<MatType>(nvirt_, nclosed_);
  copy_n(ptr_vc(), nvirt_*nclosed_, out->data());
  return out;
}


// unpack to Matrix
template<typename DataType>
shared_ptr<typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type>
  RotationMatrix<DataType>::unpack(const DataType a) const {

  const int nocc = nclosed_ + nact_;
  const int nbasis = nclosed_ + nact_ + nvirt_;
  auto out = make_shared<MatType>(nbasis, nbasis);
  fill_n(out->data(), out->size(), a);

  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc, i+nclosed_) = ele_va(j, i);
    }
    for (int j = 0; j != nclosed_; ++j) {
      out->element(i+nclosed_, j) = ele_ca(j, i);
    }
  }
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc, i) = ele_vc(j, i);
    }
  }
  for (int i = 0; i != nbasis; ++i) {
    for (int j = 0; j <= i; ++j) {
      out->element(j, i) = - detail::conj(out->element(i, j));
    }
  }
  return out;
}


template<typename DataType>
shared_ptr<typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type>
  RotationMatrix<DataType>::unpack_sym(const DataType a) const {

  const int nocc = nclosed_ + nact_;
  const int nbasis = nclosed_ + nact_ + nvirt_;
  auto out = make_shared<MatType>(nbasis, nbasis);
  fill_n(out->data(), out->size(), a);
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc, i+nclosed_) = ele_va(j, i);
    }
    for (int j = 0; j != nclosed_; ++j) {
      out->element(i+nclosed_, j) = ele_ca(j, i);
    }
  }
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc, i) = ele_vc(j, i);
    }
  }
  for (int i = 0; i != nbasis; ++i) {
    for (int j = 0; j <= i; ++j) {
      out->element(j, i) = detail::conj(out->element(i, j));
    }
  }
  return out;
}


// print matrix
template<typename DataType>
void RotationMatrix<DataType>::print(const string in) const {
  cout << "++++ " + in + " ++++" << endl;
  if (nact_ && nclosed_) {
    cout << " printing closed-active block" << endl;
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nclosed_; ++j) {
        cout << setw(10) << setprecision(6) << ele_ca(j,i);
      }
      cout << endl;
    }
  }
  if (nact_ && nvirt_) {
    cout << " printing virtual-active block" << endl;
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nvirt_; ++j) {
        cout << setw(10) << setprecision(6) << ele_va(j,i);
      }
      cout << endl;
    }
  }
  if (nclosed_ && nvirt_) {
    cout << " printing virtual-closed block" << endl;
    for (int i = 0; i != nclosed_; ++i) {
      for (int j = 0; j != nvirt_; ++j) {
        cout << setw(10) << setprecision(6) << ele_vc(j,i);
      }
      cout << endl;
    }
  }
}


template class bagel::RotationMatrix<double>;
template class bagel::RotationMatrix<std::complex<double>>;
