//
// BAGEL - Parallel electron correlation program.
// Filename: zdvec.cc
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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
#include <src/zfci/zdvec.h>

using namespace std;
using namespace bagel;

ZDvec::ZDvec(shared_ptr<const Determinants> det, const size_t ij) : det_(det), lena_(det->lena()), lenb_(det->lenb()), ij_(ij) {
  // data should be in a consecutive area to call dgemm.
  data_ = unique_ptr<complex<double>[]>(new std::complex<double>[lenb_*lena_*ij_]);
  complex<double>* tmp = data_.get();
  for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_) {
    dvec_.push_back(make_shared<ZCivec>(det_, tmp));
  }
}


ZDvec::ZDvec(const ZDvec& o) : det_(o.det_), lena_(o.lena_), lenb_(o.lenb_), ij_(o.ij_) {
  data_ = unique_ptr<complex<double>[]>(new complex<double>[lena_*lenb_*ij_]);
  complex<double>* tmp = data_.get();
  for (int i = 0; i != ij_; ++i, tmp += lenb_*lena_) {
    dvec_.push_back(make_shared<ZCivec>(det_, tmp));
  }
  copy_n(o.data(), lena_*lenb_*ij_, data());
}


ZDvec::ZDvec(shared_ptr<const ZCivec> e, const size_t ij) : det_(e->det()), lena_(e->lena()), lenb_(e->lenb()), ij_(ij) {
  data_ = unique_ptr<complex<double>[]>(new complex<double>[lena_*lenb_*ij]);
  complex<double>* tmp = data();
  for (int i = 0; i != ij; ++i, tmp += lenb_*lena_) {
    auto c = make_shared<ZCivec>(det_, tmp);
    copy_n(e->data(), lenb_*lena_, c->data());
    dvec_.push_back(c);
  }
}


// I think this is very confusiong... this is done this way in order not to delete Civec when Dvec is deleted.
ZDvec::ZDvec(shared_ptr<const ZDvec> o) : det_(o->det_), lena_(o->lena_), lenb_(o->lenb_), ij_(o->ij_) {
  for (int i = 0; i != ij_; ++i) {
    dvec_.push_back(make_shared<ZCivec>(*(o->data(i))));
  }
}


ZDvec::ZDvec(vector<shared_ptr<ZCivec>> o) : det_(o.front()->det()), ij_(o.size()) {
  lena_ = det_->lena();
  lenb_ = det_->lenb();
  for (int i = 0; i != ij_; ++i) {
    dvec_.push_back(make_shared<ZCivec>(*(o.at(i))));
  }
}


// returns a vector of Civec's which correspond to an unconverged state
vector<shared_ptr<ZCivec>> ZDvec::dvec(const vector<int>& conv) {
  vector<shared_ptr<ZCivec>> out;
  int i = 0;
  for (auto& iter : dvec_)
    if (conv[i++] == 0) out.push_back(iter);
  return out;
}


vector<shared_ptr<const ZCivec>> ZDvec::dvec(const vector<int>& conv) const {
  vector<shared_ptr<const ZCivec>> out;
  int i = 0;
  for (auto& iter : dvec_)
    if (conv[i++] == 0) out.push_back(iter);
  return out;
}


void ZDvec::set_det(shared_ptr<const Determinants> o) const {
  det_ = o;
  for (auto& i : dvec_)
    i->set_det(o);
}


complex<double> ZDvec::dot_product(const ZDvec& other) const {
  assert(ij() == other.ij());
  complex<double> sum = 0.0;
  auto j = other.dvec_.begin();
  for (auto& i : dvec_) {
    sum += i->dot_product(**j);
    ++j;
  }
  return sum;
}


void ZDvec::ax_plus_y(complex<double> a, const ZDvec& other) {
  auto j = other.dvec_.begin();
  for (auto& i : dvec_) {
    i->ax_plus_y(a, **j);
    ++j;
  }
}


void ZDvec::scale(const complex<double> a) {
  for (auto& i : dvec_)
    i->scale(a);
}


ZDvec& ZDvec::operator/=(const ZDvec& o) {
  assert(dvec().size() == o.dvec().size());
  auto j = o.dvec().begin();
  for (auto i = dvec().begin(); i != dvec().end(); ++i, ++j)
    **i /= **j;
  return *this;
}


ZDvec ZDvec::operator/(const ZDvec& o) const {
  ZDvec out(*this);
  out /= o;
  return out;
}


void ZDvec::orthog(shared_ptr<const ZDvec> o) {
  if (o->ij() != ij()) throw logic_error("ZDvec::orthog called inconsistently");
  auto j = o->dvec().begin();
  for (auto i = dvec().begin(); i != dvec().end(); ++i, ++j)
    (*i)->orthog(*j);
}


void ZDvec::project_out(shared_ptr<const ZDvec> o) {
  if (o->ij() != ij()) throw logic_error("ZDvec::project_out called inconsistently");
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

shared_ptr<ZDvec> ZDvec::spin() const {
  vector<shared_ptr<ZCivec>> ccvec;
  for (auto& cc : dvec_) {
    ccvec.push_back(cc->spin());
  }

  return make_shared<ZDvec>(ccvec);
}

shared_ptr<ZDvec> ZDvec::spinflip(shared_ptr<const Determinants> det) const {
  if(det == nullptr) det = det_->transpose();

  vector<shared_ptr<ZCivec>> ccvec;
  for (auto& cc : dvec_) {
    ccvec.push_back(cc->transpose(det));
  }

  return make_shared<ZDvec>(ccvec);
}

shared_ptr<ZDvec> ZDvec::spin_lower(shared_ptr<const Determinants> det) const {
  if (det == nullptr)
    det = make_shared<Determinants>(det_->norb(), det_->nelea()-1, det_->neleb()+1, det_->compress(), true);

  vector<shared_ptr<ZCivec>> ccvec;
  for (auto& cc : dvec_) {
    ccvec.push_back(cc->spin_lower(det));
  }

  return make_shared<ZDvec>(ccvec);
}

shared_ptr<ZDvec> ZDvec::spin_raise(shared_ptr<const Determinants> det) const {
  if(det == nullptr)
    det = make_shared<Determinants>(det_->norb(), det_->nelea()+1, det_->neleb()-1, det_->compress(), true);

  vector<shared_ptr<ZCivec>> ccvec;
  for (auto& cc : dvec_) {
    ccvec.push_back(cc->spin_raise(det));
  }

  return make_shared<ZDvec>(ccvec);
}

void ZDvec::print(const double thresh) const {
  int j = 0;
  for (auto& iter : dvec_) {
    cout << endl;
    cout << "     * ci vector, state " << setw(3) << j++ << endl;
    // TODO spin_expectation is not correct for relativistic cases << ", <S^2> = " << setw(6) << setprecision(4) << iter->spin_expectation() << endl;
    iter->print(thresh);
  }
}
