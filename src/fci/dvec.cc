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
#include <src/fci/dvec.h>

using namespace std;
using namespace bagel;


template<>
shared_ptr<Dvector<double>> Dvector<double>::spin() const {
  vector<shared_ptr<Civec>> ccvec;
  for (auto& cc : dvec_) {
    ccvec.push_back(cc->spin());
  }

  return make_shared<Dvector<double>>(ccvec);
}

template<>
shared_ptr<Dvector<double>> Dvector<double>::spinflip(shared_ptr<const Determinants> det) const {
  if(det == nullptr) det = det_->transpose();

  vector<shared_ptr<Civec>> ccvec;
  for (auto& cc : dvec_) {
    ccvec.push_back(cc->transpose(det));
  }

  return make_shared<Dvector<double>>(ccvec);
}

template<>
shared_ptr<Dvector<double>> Dvector<double>::spin_lower(shared_ptr<const Determinants> det) const {
  if (det == nullptr)
    det = make_shared<Determinants>(det_->norb(), det_->nelea()-1, det_->neleb()+1, det_->compress(), true);

  vector<shared_ptr<Civec>> ccvec;
  for (auto& cc : dvec_) {
    ccvec.push_back(cc->spin_lower(det));
  }

  return make_shared<Dvector<double>>(ccvec);
}

template<>
shared_ptr<Dvector<double>> Dvector<double>::spin_raise(shared_ptr<const Determinants> det) const {
  if(det == nullptr)
    det = make_shared<Determinants>(det_->norb(), det_->nelea()+1, det_->neleb()-1, det_->compress(), true);

  vector<shared_ptr<Civec>> ccvec;
  for (auto& cc : dvec_) {
    ccvec.push_back(cc->spin_raise(det));
  }

  return make_shared<Dvector<double>>(ccvec);
}

