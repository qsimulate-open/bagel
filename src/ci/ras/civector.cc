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

#include <src/util/taskqueue.h>
#include <src/ci/ras/civector.h>
#include <src/ci/ras/civec_spinop.h>
#include <src/ci/ras/apply_block.h>

using namespace std;
using namespace bagel;


template<typename DataType>
RASCivector<DataType>::RASCivector(shared_ptr<const RASDeterminants> det) : RASCivector_impl<DataType>(det) {
  data_ = unique_ptr<DataType[]>(new DataType[size()]);
  fill_n(data_.get(), size(), 0.0);
  size_t sz = 0;
  for (auto& ipair : det->blockinfo())
    if (!ipair->empty()) {
      blocks_.push_back(make_shared<RBlock>(ipair->stringsa(), ipair->stringsb(), data_.get()+sz, sz));
      sz += blocks_.back()->size();
    } else {
      blocks_.push_back(nullptr);
    }
}


template<typename DataType>
RASCivector<DataType>& RASCivector<DataType>::operator=(RASCivector<DataType>&& o) {
  assert(*o.det() == *det());
  data_ = move(o.data_);
  blocks_ = move(o.blocks_);
  return *this;
}


template<typename DataType>
shared_ptr<RASCivector<DataType>> RASCivector<DataType>::apply(const int orbital, const bool action, const bool spin) const {
  // action: true -> create; false -> annihilate
  // spin: true -> alpha; false -> beta
  shared_ptr<const RASDeterminants> sdet = this->det();

  const int ras1 = sdet->ras(0);
  const int ras2 = sdet->ras(1);
  const int ras3 = sdet->ras(2);

  // 0 -> RASI, 1 -> RASII, 2 -> RASIII
  const int ras_space = ( orbital >= ras1 ) + (orbital >= ras1 + ras2);

  auto to_array = [] (shared_ptr<const RASBlock<DataType>> block) {
    auto sa = block->stringsa();
    auto sb = block->stringsb();
    return array<int, 6>({sa->nholes(), sb->nholes(), sa->nele2(), sb->nele2(), sa->nparticles(), sb->nparticles()});
  };

  auto op_on_array = [&ras_space, &action, &spin] ( array<int, 6> in ) {
    const int mod = ( action ? +1 : -1 ) * ( ras_space == 0 ? -1 : 1 );
    array<int, 6> out = in;
    out[2*ras_space] += (spin ? mod : 0);
    out[2*ras_space+1] += (spin ? 0 : mod);
    return out;
  };

  RAS::Apply_block apply_block(orbital, action, spin);

  const int mod = action ? +1 : -1;
  const int telea = sdet->nelea() + (spin ? mod : 0);
  const int teleb = sdet->neleb() + (spin ? 0 : mod);
  const int tholes = max(sdet->max_holes() - ((ras_space == 0) ? mod : 0), 0);
  const int tparts = max(sdet->max_particles() + ((ras_space == 2) ? mod : 0), 0);

  auto tdet = make_shared<const RASDeterminants>(ras1, ras2, ras3, telea, teleb, tholes, tparts, true);
  auto out = make_shared<RASCivector<DataType>>(tdet);

  for (shared_ptr<const RASBlock<double>> soblock : this->blocks()) {
    if (!soblock) continue;
    array<int, 6> tar_array = op_on_array(to_array(soblock));
    if (all_of(tar_array.begin(), tar_array.end(), [] (int i) { return i >= 0; })) {
      shared_ptr<RASBlock<double>> tarblock = out->block(tar_array[0], tar_array[1], tar_array[4], tar_array[5]);
      if (tarblock) apply_block(soblock, tarblock, false);
    }
  }
  return out;
}


template<typename DataType>
RASCivecView_<DataType>::RASCivecView_(shared_ptr<const RASDeterminants> det, double* const data)
 : RASCivector_impl<DataType>(det), data_ptr_(data), can_write_(true) {
  size_t sz = 0;
  for (auto& ipair : det->blockinfo()) {
    if (!ipair->empty()) {
      blocks_.push_back(make_shared<RBlock>(ipair->stringsa(), ipair->stringsb(), data+sz, sz));
      sz += blocks_.back()->size();
    } else {
      blocks_.push_back(nullptr);
    }
  }
}


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
  assert((tdet->nelea() == sdet->nelea()-1) && (tdet->neleb() == sdet->neleb()+1));
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_lower_impl(*this, *out);
  return out;
}

template<> shared_ptr<RASCivector<double>> RASCivecView_<double>::spin_lower(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()-1, sdet->neleb()+1);
  assert((tdet->nelea() == sdet->nelea()-1) && (tdet->neleb() == sdet->neleb()+1));
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_lower_impl(*this, *out);
  return out;
}


template<> shared_ptr<RASCivector<double>> RASCivector<double>::spin_raise(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()+1, sdet->neleb()-1);
  assert((tdet->nelea() == sdet->nelea()+1) && (tdet->neleb() == sdet->neleb()-1));
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_raise_impl(*this, *out);
  return out;
}

template<> shared_ptr<RASCivector<double>> RASCivecView_<double>::spin_raise(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()+1, sdet->neleb()-1);
  assert((tdet->nelea() == sdet->nelea()+1) && (tdet->neleb() == sdet->neleb()-1));
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_raise_impl(*this, *out);
  return out;
}


template class bagel::RASCivector<double>;
template class bagel::RASCivecView_<double>;
