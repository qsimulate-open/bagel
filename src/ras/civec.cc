//
// BAGEL - Parallel electron correlation program.
// Filename: ras/civector.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#include <iomanip>
#include <unordered_map>

#include <src/ras/civector.h>
#include <src/util/taskqueue.h>
#include <src/ras/civec_spinop.h>

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
