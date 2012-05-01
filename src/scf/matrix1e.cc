//
// Newint - Parallel electron correlation program.
// Filename: matrix1e.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/scf/matrix1e.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <cassert>
#include <cmath>

using namespace std;

typedef std::shared_ptr<Geometry> RefGeometry;
typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<Shell> RefShell;

Matrix1e::Matrix1e(const RefGeometry geom) : data_(new double[geom->nbasis()*geom->nbasis()]), geom_(geom), nbasis_(geom->nbasis()) {
  mdim_ = ndim_ = nbasis_;
  zero();
}


Matrix1e::Matrix1e(const RefGeometry geom, const int n, const int m)
 : data_(new double[geom->nbasis()*geom->nbasis()]), geom_(geom), nbasis_(geom->nbasis()) {
  ndim_ = n;
  mdim_ = m;
  zero();
}


Matrix1e::Matrix1e(const Matrix1e& o)
 : data_(new double[o.geom_->nbasis()*o.geom_->nbasis()]), geom_(o.geom_), nbasis_(o.nbasis_), ndim_(o.ndim_), mdim_(o.mdim_) {
  copy(o.data(), o.data() + nbasis_*nbasis_, data());
}


void Matrix1e::init() {

  const vector<RefAtom> atoms = geom_->atoms(); 
  vector<Atom>::const_iterator aiter0, aiter1;
  vector<RefShell>::const_iterator biter0, biter1;

  const vector<vector<int> > offsets = geom_->offsets();

  const int nbasis = geom_->nbasis();

  // only lower half will be stored
  for (int iatom0 = 0; iatom0 != geom_->natom(); ++iatom0) {
    // iatom1 = iatom1;
    const RefAtom catom0 = atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = offsets[iatom0];
    const vector<RefShell> shell0 = catom0->shells();

    for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
      const int offset0 = coffset0[ibatch0]; 
      RefShell b0 = shell0[ibatch0];
      for (int ibatch1 = ibatch0; ibatch1 != numshell0; ++ibatch1) {
        const int offset1 = coffset0[ibatch1]; 
        RefShell b1 = shell0[ibatch1];
        vector<RefShell> input;
        input.push_back(b1);
        input.push_back(b0);

        computebatch(input, offset0, offset1, nbasis);
      }
    } 

    for (int iatom1 = iatom0 + 1; iatom1 != geom_->natom(); ++iatom1) {
      const RefAtom catom1 = atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = offsets[iatom1];
      const vector<RefShell> shell1 = catom1->shells();

      for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
        const int offset0 = coffset0[ibatch0]; 
        RefShell b0 = shell0[ibatch0];
        for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
          const int offset1 = coffset1[ibatch1]; 
          RefShell b1 = shell1[ibatch1];
          vector<RefShell> input;
          input.push_back(b1);
          input.push_back(b0);

          computebatch(input, offset0, offset1, nbasis);

        }
      } 
    }
  } 

}


Matrix1e::~Matrix1e() {
}


void Matrix1e::computebatch(const std::vector<RefShell>&, const int, const int, const int) {
  assert(false);
}


shared_ptr<Matrix1e> Matrix1e::resize(shared_ptr<Geometry> g, const int n) const {
  shared_ptr<Matrix1e> out(new Matrix1e(g));
  out->ndim_ = n; 
  out->mdim_ = mdim_;
  assert(n <= g->nbasis());
  for (int i = 0; i != mdim_; ++i) {
    for (int j = 0; j != ndim_; ++j) {
      out->data_[j+i*g->nbasis()] = data_[j+i*nbasis_];
    }
  }
  return out;
}


shared_ptr<Matrix1e> Matrix1e::slice(const int start, const int fence) const {
  shared_ptr<Matrix1e> out(new Matrix1e(geom_));
  out->ndim_ = ndim_;
  out->mdim_ = fence - start;
  assert(fence <= geom_->nbasis());

  copy(data_.get()+start*ndim_, data_.get()+fence*ndim_, out->data_.get());
  return out;
}

shared_ptr<Matrix1e> Matrix1e::merge(const shared_ptr<const Matrix1e> o) const {
  shared_ptr<Matrix1e> out(new Matrix1e(*this));
  assert(nbasis_ == o->geom()->nbasis());
  assert(ndim_ == o->ndim_);
  out->ndim_ = ndim_;
  out->mdim_ = mdim_ + o->mdim_;
  assert(out->mdim_ <= nbasis_);

  copy(o->data_.get(), o->data_.get()+nbasis_*o->mdim_, out->data_.get()+mdim_*nbasis_);
  return out;
}


void Matrix1e::symmetrize() {
  for (int i = 0; i != nbasis_; ++i) {
    for (int j = i + 1; j != nbasis_; ++j) {
      data_[i + j * nbasis_] = data_[j + i * nbasis_]; 
    }
  } 
}


Matrix1e Matrix1e::operator+(const Matrix1e& o) const {
  Matrix1e out(geom_);
  out.ndim_ = ndim_;
  out.mdim_ = mdim_;
  const int size = nbasis_ * nbasis_;
  const double* odata = o.data();
  double* outdata = out.data();

  daxpy_(size, 1.0, data(), 1, outdata, 1); 
  daxpy_(size, 1.0, odata, 1, outdata, 1); 

  return out; 
}


Matrix1e& Matrix1e::operator+=(const Matrix1e& o) {
  daxpy_(nbasis_*nbasis_, 1.0, o.data(), 1, data(), 1); 
  return *this; 
}


Matrix1e& Matrix1e::operator-=(const Matrix1e& o) {
  daxpy_(nbasis_*nbasis_, -1.0, o.data(), 1, data(), 1); 
  return *this; 
}


Matrix1e& Matrix1e::operator=(const Matrix1e& o) {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  dcopy_(nbasis_*nbasis_, o.data(), 1, data(), 1); 
  return *this; 
}


Matrix1e Matrix1e::operator-(const Matrix1e& o) const {
  Matrix1e out(geom_);
  out.ndim_ = ndim_;
  out.mdim_ = mdim_;
  const int size = nbasis_ * nbasis_;
  const double* odata = o.data();
  double* outdata = out.data();
  daxpy_(size, 1.0, data(), 1, outdata, 1); 
  daxpy_(size, -1.0, odata, 1, outdata, 1); 

  return out; 
}


Matrix1e Matrix1e::operator*(const Matrix1e& o) const {
  Matrix1e out(geom_);
  const int l = ndim_;
  const int m = mdim_; 
  assert(mdim_ == o.ndim());
  const int n = o.mdim(); 
  const double* odata = o.data();
  double* outdata = out.data();

  dgemm_("N", "N", l, n, m, 1.0, data(), nbasis_, odata, nbasis_, 0.0, outdata, nbasis_); 

  out.ndim_ = l;
  out.mdim_ = n; 
  return out;
}


Matrix1e& Matrix1e::operator*=(const Matrix1e& o) {
  Matrix1e out(geom_);
  const int l = ndim_;
  const int m = mdim_; 
  assert(mdim_ == o.ndim());
  const int n = o.mdim(); 
  const double* odata = o.data();
  dgemm_("N", "N", l, n, m, 1.0, data(), nbasis_, odata, nbasis_, 0.0, out.data(), nbasis_); 
  *this = out;
}


Matrix1e Matrix1e::operator*(const double& a) const {
  Matrix1e out(*this);
  dscal_(nbasis_*nbasis_, a, out.data(), 1);
  return out;
}
Matrix1e Matrix1e::operator/(const double& a) const {
  Matrix1e out(*this);
  dscal_(nbasis_*nbasis_, 1.0/a, out.data(), 1);
  return out;
}

Matrix1e& Matrix1e::operator*=(const double& a) {
  dscal_(nbasis_*nbasis_, a, data_, 1);
  return *this;
}
Matrix1e& Matrix1e::operator/=(const double& a) {
  dscal_(nbasis_*nbasis_, 1.0/a, data_, 1);
  return *this;
}


Matrix1e Matrix1e::operator%(const Matrix1e& o) const {
  Matrix1e out(geom_);
  const int l = mdim_;
  const int m = ndim_; 
  assert(ndim_ == o.ndim());
  const int n = o.mdim(); 
  const double* odata = o.data();
  double* outdata = out.data();

  dgemm_("T", "N", l, n, m, 1.0, data(), nbasis_, odata, nbasis_, 0.0, outdata, nbasis_); 

  out.ndim_ = l;
  out.mdim_ = n;
    
  return out;
}


Matrix1e Matrix1e::operator^(const Matrix1e& o) const {
  Matrix1e out(geom_);
  const int l = ndim_;
  const int m = mdim_; 
  assert(mdim_ == o.mdim());
  const int n = o.ndim(); 
  const double* odata = o.data();
  double* outdata = out.data();

  dgemm_("N", "T", l, n, m, 1.0, data(), nbasis_, odata, nbasis_, 0.0, outdata, nbasis_); 

  out.ndim_ = l;
  out.mdim_ = n;
    
  return out;
}


void Matrix1e::diagonalize(double* eig) {
  // assume that the matrix is symmetric
  // the leading order (nbasis supplied)
  
  int info;
  const int lwork = nbasis_*6;
  unique_ptr<double[]> work(new double[lwork]);
  dsyev_("V", "L", ndim_, data(), nbasis_, eig, work.get(), lwork, info); 

  if(info) throw runtime_error("diagonalize failed");

}


void Matrix1e::svd(shared_ptr<Matrix1e> U, shared_ptr<Matrix1e> V) {
  assert(U->ndim() == ndim_ && U->mdim() == ndim_);
  assert(V->ndim() == mdim_ && V->mdim() == mdim_);
  const int lwork = 10*max(ndim_, mdim_);
  unique_ptr<double[]> work(new double[lwork]);
  unique_ptr<double[]> S(new double[min(ndim_, mdim_)]);
/*
  SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
 $                   WORK, LWORK, INFO )
 */
  double* cblock = data();
  double* ublock = U->data();
  double* vblock = V->data();
  int info = 0;
  dgesvd_("A", "A", &ndim_, &mdim_, cblock, &ndim_, S.get(), ublock, &ndim_, vblock, &mdim_, work.get(), &lwork, &info);
  if (info != 0) throw runtime_error("dgesvd failed in Matrix1e::svd");
}


void Matrix1e::daxpy(const double a, const Matrix1e& o) {
  daxpy_(nbasis_*nbasis_, a, o.data(), 1, data(), 1); 
}


void Matrix1e::daxpy(const double a, const std::shared_ptr<const Matrix1e> o) {
  daxpy_(nbasis_*nbasis_, a, o->data(), 1, data(), 1); 
}


const double Matrix1e::ddot(const Matrix1e& o) const {
  return ddot_(nbasis_*nbasis_, data(), 1, o.data(), 1); 
}


const double Matrix1e::ddot(const std::shared_ptr<const Matrix1e> o) const {
  return ddot_(nbasis_*nbasis_, data(), 1, o->data(), 1); 
}


const double Matrix1e::rms() const {
  return ::sqrt(ddot(*this) / (ndim_ * mdim_));
}


const double Matrix1e::trace() const {
  double out = 0.0;
  for (int i = 0; i != ndim_; ++i) out += data_[i * nbasis_ + i]; 
  return out;
}


shared_ptr<Matrix1e> Matrix1e::exp(const int deg) const {
  shared_ptr<Matrix1e> out(new Matrix1e(geom_, ndim_, mdim_));
  Matrix1e buf(*this);
  assert(ndim_ == mdim_);

  for (int i = deg; i != 1; --i) {
    const double inv = 1.0/static_cast<double>(i);
    buf *= inv;
    for (int j = 0; j != ndim_; ++j) buf.element(j,j) += 1.0; 
    *out = (*this)*buf;
    if (i != 1) buf = *out;
  }
  for (int j = 0; j != ndim_; ++j) out->element(j,j) += 1.0; 
  return out;
}


shared_ptr<Matrix1e> Matrix1e::log(const int deg) const {
  shared_ptr<Matrix1e> out(new Matrix1e(geom_, ndim_, mdim_));
  Matrix1e buf(*this);
  for (int j = 0; j != ndim_; ++j) buf.element(j,j) -= 1.0;
  assert(ndim_ == mdim_);

  for (int i = deg; i != 1; --i) {
    const double inv = -static_cast<double>(i-1)/static_cast<double>(i);
    buf *= inv;
    for (int j = 0; j != ndim_; ++j) buf.element(j,j) += 1.0; 
    *out = (*this)*buf - buf;
    if (i != 1) buf = *out;
  }
  return out;
}


unique_ptr<double[]> Matrix1e::diag() const {
  if (ndim_ != mdim_) throw logic_error("illegal call of Matrix1e::diag()");
  unique_ptr<double[]> out(new double[ndim_]);
  for (int i = 0; i != ndim_; ++i) {
    out[i] = element(i,i); 
  }
  return move(out);
}


shared_ptr<Matrix1e> Matrix1e::transpose() const {
  shared_ptr<Matrix1e> out(new Matrix1e(*this));
  mytranspose_(data(), &nbasis_, &nbasis_, out->data());
  return out;
}


void Matrix1e::purify_unitary() {
  assert(ndim_ == mdim_);
#if 0
  Matrix1e buf(*this ^ *this);
  const int lwork = 5*ndim_;
  unique_ptr<double[]> work(new double[lwork]);
  unique_ptr<double[]> vec(new double[ndim_]);
  int info;
  dsyev_("V", "U", ndim_, buf.data(), ndim_, vec.get(), work.get(), lwork, info); 
  if (info) throw runtime_error("dsyev failed in Matrix1e::purify_unitary");
  if (vec[0] < 0.95)        cout << "   --- smallest eigenvalue in purify_unitary() " << vec[0] << endl;
  if (vec[ndim_-1] > 1.05)  cout << "   --- largest eigenvalue in purify_unitary() " << vec[ndim_-1] << endl;
  for (int i = 0; i != ndim_; ++i) {
    for (int j = 0; j != ndim_; ++j) {
      buf.element(j,i) /= std::sqrt(std::sqrt(vec[i]));
    }
  }
  *this = ((buf ^ buf) * *this);
#else
  // Schmidt orthogonalization
  for (int i = 0; i != ndim_; ++i) {
    for (int j = 0; j != i; ++j) {
      const double a = ddot_(ndim_, &data_[i*ndim_], 1, &data_[j*ndim_], 1);
      daxpy_(ndim_, -a, &data_[j*ndim_], 1, &data_[i*ndim_], 1);
    }
    const double b = 1.0/sqrt(ddot_(ndim_, &data_[i*ndim_], 1, &data_[i*ndim_], 1));
    dscal_(ndim_, b, &data_[i*ndim_], 1);
  }
#endif

}

void Matrix1e::purify_redrotation(const int nclosed, const int nact, const int nvirt) {

#if 1
  for (int g = 0; g != nclosed; ++g)
    for (int h = 0; h != nclosed; ++h)
      element(h,g)=0.0;
  for (int g = 0; g != nact; ++g)
    for (int h = 0; h != nact; ++h)
      element(h+nclosed,g+nclosed)=0.0;
  for (int g = 0; g != nvirt; ++g)
    for (int h = 0; h != nvirt; ++h)
      element(h+nclosed+nact,g+nclosed+nact)=0.0;
  for (int i = 0; i != nbasis_; ++i) {
    for (int j = 0; j != i; ++j) {
      const double ele = (element(j,i) - element(i,j)) * 0.5;
      element(j,i) = ele; 
      element(i,j) = -ele; 
    }
  }
#endif

}

void Matrix1e::purify_idempotent(const Matrix1e& s) {
  *this = *this * s * *this * 3.0 - *this * s * *this * s * *this * 2.0; 
}


// in-place matrix inverse (practically we use buffer area)
void Matrix1e::inverse() {
  assert(ndim_ == mdim_);
  shared_ptr<Matrix1e> buf = this->clone();
  buf->unit();

  const int lwork = 3*ndim_;
  int info;
  unique_ptr<int[]> ipiv(new int[ndim_]);
  dgesv_(&ndim_, &ndim_, data(), &ndim_, ipiv.get(), buf->data(), &ndim_, &info);
  if (info) throw runtime_error("dsysv failed in Matrix1e::inverse()");

  copy(buf->data(), buf->data()+nbasis_*nbasis_, data());

}


double Matrix1e::orthog(const std::list<std::shared_ptr<const Matrix1e> > o) {
  for (auto iter = o.begin(); iter != o.end(); ++iter) {
    const double m = this->ddot(*iter);
    this->daxpy(-m, *iter);
  }
  const double n = norm();
  *this /= n;
  return n;
}

void Matrix1e::print(const string name, const int size) const { 
 
  cout << "++++ " + name + " ++++" << endl;
  for (int i = 0; i != min(size,nbasis_); ++i) {
    for (int j = 0; j != min(size,nbasis_); ++j) {
      cout << fixed << setw(9) << setprecision(6) << data_[j * nbasis_ + i]  << " "; 
    }
    cout << endl;
  }

} 

