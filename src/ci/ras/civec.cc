//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/civector.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#include <iomanip>
#include <unordered_map>

#include <src/ci/ras/civector.h>
#include <src/util/taskqueue.h>
#include <src/ci/ras/civec_spinop.h>

template class bagel::RASCivector<double>;
template class bagel::RASCivecView_<double>;

using namespace std;
using namespace bagel;

template<>
shared_ptr<RASCivector<double>> RASCivector<double>::spin() const {
  auto out = make_shared<RASCivector<double>>(det_);
  RAS::spin_impl(*this, *out);
  return out;
}

template<>
shared_ptr<RASCivector<double>> RASCivecView_<double>::spin() const {
  auto out = make_shared<RASCivector<double>>(det_);
  RAS::spin_impl(*this, *out);
  return out;
}

template<> shared_ptr<RASCivector<double>> RASCivector<double>::spin_lower(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()-1, sdet->neleb()+1);
  assert( (tdet->nelea() == sdet->nelea()-1) && (tdet->neleb() == sdet->neleb()+1) );
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_lower_impl(*this, *out);
  return out;
}

template<> shared_ptr<RASCivector<double>> RASCivecView_<double>::spin_lower(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()-1, sdet->neleb()+1);
  assert( (tdet->nelea() == sdet->nelea()-1) && (tdet->neleb() == sdet->neleb()+1) );
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_lower_impl(*this, *out);
  return out;
}

template<> shared_ptr<RASCivector<double>> RASCivector<double>::spin_raise(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()+1, sdet->neleb()-1);
  assert( (tdet->nelea() == sdet->nelea()+1) && (tdet->neleb() == sdet->neleb()-1) );
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_raise_impl(*this, *out);
  return out;
}

template<> shared_ptr<RASCivector<double>> RASCivecView_<double>::spin_raise(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()+1, sdet->neleb()-1);
  assert( (tdet->nelea() == sdet->nelea()+1) && (tdet->neleb() == sdet->neleb()-1) );
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_raise_impl(*this, *out);
  return out;
}
