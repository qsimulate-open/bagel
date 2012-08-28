//
// BAGEL - Parallel electron correlation program.
// Filename: civec.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/newfci/civec.h>

using namespace std;

NewCivec::NewCivec(shared_ptr<const NewDeterminants> det) : det_(det), lena_(det->lena()), lenb_(det->lenb()) {
  cc_ = unique_ptr<double[]>(new double[lena_*lenb_]);
  cc_ptr_ = cc_.get();
  fill(cc(), cc() + lena_*lenb_, 0.0);
}


NewCivec::NewCivec(shared_ptr<const NewDeterminants> det, double* din_) : det_(det), lena_(det->lena()), lenb_(det->lenb()) {
  cc_ = unique_ptr<double[]>(new double[lena_*lenb_]);
  cc_ptr_ = din_;
  fill(cc(), cc() + lena_*lenb_, 0.0);
}


NewCivec::NewCivec(const NewCivec& o) : det_(o.det_), lena_(o.lena_), lenb_(o.lenb_) {
  cc_ = unique_ptr<double[]>(new double[lena_*lenb_]);
  cc_ptr_ = cc_.get();
  copy(o.cc(), o.cc() + lena_*lenb_, cc());
}


NewCivec::NewCivec(shared_ptr<NewCivec> o, shared_ptr<const NewDeterminants> det) : det_(det), lena_(o->lena_), lenb_(o->lenb_) {
  assert(lena_ == det->lena() && lenb_ == det->lenb());
  cc_ = move(o->cc_);
  cc_ptr_ = cc_.get();
}


shared_ptr<NewCivec> NewCivec::transpose() const {
  shared_ptr<NewCivec> ct(new NewCivec(det_));
  double* cct = ct->data(); 
  mytranspose_(cc(), &lenb_, &lena_, cct); 
  return ct;
}


double NewCivec::ddot(const NewCivec& other) const {
  return ddot_(lena_*lenb_, cc(), 1, other.data(), 1);
}


void NewCivec::daxpy(double a, const NewCivec& other) {
  daxpy_(lena_*lenb_, a, other.data(), 1, cc(), 1);
}


double NewCivec::norm() const {
  return sqrt(ddot_(lena_*lenb_, cc(), 1, cc(), 1));
}


void NewCivec::scale(const double a) {
  dscal_(lena_*lenb_, a, cc(), 1);
}


double NewCivec::variance() const {
  return ddot_(lena_*lenb_, cc(), 1, cc(), 1) / (lena_*lenb_);
}


double NewCivec::orthog(list<shared_ptr<const NewCivec> > c) {
  for (auto iter = c.begin(); iter != c.end(); ++iter)
    project_out(*iter);
  const double norm = this->norm();
  const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
  dscal_(lena_*lenb_, scal, cc(), 1);
  return 1.0/scal; 
}

double NewCivec::orthog(shared_ptr<const NewCivec> o) {
  list<shared_ptr<const NewCivec> > v; v.push_back(o);
  return orthog(v);
}


NewCivec& NewCivec::operator/=(const NewCivec& o) {
  for (int i = 0; i != size(); ++i) {
    data(i) /= o.data(i);
  }
  return *this;
}

NewCivec NewCivec::operator/(const NewCivec& o) const {
  NewCivec out(*this);
  out /= o;
  return out;
}
