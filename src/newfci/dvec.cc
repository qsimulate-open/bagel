//
// BAGEL - Parallel electron correlation program.
// Filename: dvec.cc
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

#include <stdexcept>
#include <src/newfci/dvec.h>

using namespace std;

NewDvec::NewDvec(shared_ptr<const NewDeterminants> det, const size_t ij) : det_(det), lena_(det->lena()), lenb_(det->lenb()), ij_(ij) {
  // data should be in a consecutive area to call dgemm.
  data_ = unique_ptr<double[]>(new double[lenb_*lena_*ij_]);
  double* tmp = data_.get();
  for (int i = 0; i != ij_; ++i, tmp+=lenb_*lena_) {
    shared_ptr<NewCivec> c(new NewCivec(det_, tmp)); 
    dvec_.push_back(c);
  }
}


NewDvec::NewDvec(const NewDvec& o) : det_(o.det_), lena_(o.lena_), lenb_(o.lenb_), ij_(o.ij_) {
  data_ = unique_ptr<double[]>(new double[lena_*lenb_*ij_]);
  double* tmp = data_.get();
  for (int i = 0; i != ij_; ++i, tmp+=lenb_*lena_) {
    shared_ptr<NewCivec> c(new NewCivec(det_, tmp)); 
    dvec_.push_back(c);
  }
  std::copy(o.data(), o.data()+lena_*lenb_*ij_, data());
}


NewDvec::NewDvec(shared_ptr<const NewCivec> e, const size_t ij) : det_(e->det()), lena_(e->lena()), lenb_(e->lenb()), ij_(ij) {
  data_ = unique_ptr<double[]>(new double[lena_*lenb_*ij]);
  double* tmp = data();
  for (int i = 0; i != ij; ++i, tmp+=lenb_*lena_) {
    shared_ptr<NewCivec> c(new NewCivec(det_, tmp)); 
    dcopy_(lenb_*lena_, e->data(), 1, c->data(), 1);
    dvec_.push_back(c);
  }
}


// I think this is very confusiong... this is done this way in order not to delete NewCivec when NewDvec is deleted. 
NewDvec::NewDvec(shared_ptr<const NewDvec> o) : det_(o->det_), lena_(o->lena_), lenb_(o->lenb_), ij_(o->ij_) {
  for (int i = 0; i != ij_; ++i) {
    shared_ptr<NewCivec> c(new NewCivec(*(o->data(i))));
    dvec_.push_back(c);
  }
}


NewDvec::NewDvec(vector<shared_ptr<NewCivec> > o) : det_(o.front()->det()), ij_(o.size()) {
  lena_ = det_->lena();
  lenb_ = det_->lenb();
  dvec_ = o;
}


// returns a vector of NewCivec's which correspond to an unconverged state
vector<shared_ptr<NewCivec> > NewDvec::dvec(const vector<int>& conv) {
  vector<shared_ptr<NewCivec> > out;
  int i = 0;
  for (auto iter = dvec_.begin(); iter != dvec_.end(); ++iter, ++i) {
    if (conv[i] == 0) out.push_back(*iter);
  }
  return out;
}


vector<shared_ptr<const NewCivec> > NewDvec::dvec(const vector<int>& conv) const {
  vector<shared_ptr<const NewCivec> > out;
  int i = 0;
  for (auto iter = dvec_.begin(); iter != dvec_.end(); ++iter, ++i) {
    if (conv[i] == 0) out.push_back(*iter);
  }
  return out;
}


void NewDvec::set_det(shared_ptr<const NewDeterminants> o) const {
  det_ = o;
  for (auto i = dvec_.begin(); i != dvec_.end(); ++i) (*i)->set_det(o); 
}


double NewDvec::ddot(const NewDvec& other) const {
  assert(ij() == other.ij());
  double sum = 0.0;
  for (auto i = dvec_.begin(), j = other.dvec_.begin(); i != dvec_.end(); ++i, ++j)
    sum += (*i)->ddot(**j); 
  return sum;
}


void NewDvec::daxpy(double a, const NewDvec& other) {
  auto j = other.dvec_.begin();
  for (auto i = dvec_.begin(); i != dvec_.end(); ++i, ++j)
    daxpy_(lena_*lenb_, a, (*j)->data(), 1, (*i)->data(), 1);
}


void NewDvec::scale(const double a) {
  for (auto i = dvec_.begin(); i != dvec_.end(); ++i)
    dscal_(lena_*lenb_, a, (*i)->data(), 1);
}


NewDvec& NewDvec::operator/=(const NewDvec& o) {
  assert(dvec().size() == o.dvec().size());
  auto j = o.dvec().begin();
  for (auto i = dvec().begin(); i != dvec().end(); ++i, ++j)
    **i /= **j; 
  return *this;
}


NewDvec NewDvec::operator/(const NewDvec& o) const {
  NewDvec out(*this);
  out /= o;
  return out;
}


void NewDvec::orthog(shared_ptr<const NewDvec> o) {
  if (o->ij() != ij()) throw logic_error("NewDvec::orthog called inconsistently"); 
  auto j = o->dvec().begin();
  for (auto i = dvec().begin(); i != dvec().end(); ++i, ++j)
    (*i)->orthog(*j);
}


void NewDvec::project_out(shared_ptr<const NewDvec> o) {
  if (o->ij() != ij()) throw logic_error("NewDvec::project_out called inconsistently"); 
#if 1
  auto j = o->dvec().begin();
  // simply project out each CI vector
  for (auto i = dvec().begin(); i != dvec().end(); ++i, ++j)
    (*i)->project_out(*j);
#else
  for (auto i = dvec().begin(); i != dvec().end(); ++i)
    for (auto j = o->dvec().begin(); j != o->dvec().end(); ++j)
      (*i)->project_out(*j);
#endif
}


void NewDvec::print(const double thresh) const {
  int j = 0;
  for (auto iter = dvec().begin(); iter != dvec().end(); ++iter, ++j) {
    cout << endl;
    cout << "     * ci vector, state " << setw(3) << j << endl; 
    (*iter)->print(thresh);
  }
}
