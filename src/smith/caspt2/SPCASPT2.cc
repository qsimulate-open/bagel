//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: SPCASPT2.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/ci/fci/fci.h>
#include <src/smith/smith_util.h>
#include <src/smith/caspt2/SPCASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

// TODO assuming a single-state calculation

SPCASPT2::SPCASPT2::SPCASPT2(const CASPT2::CASPT2& cas) {
  virt_    = cas.virt_;
  active_  = cas.active_;
  closed_  = cas.closed_;
  rvirt_   = cas.rvirt_;
  ractive_ = cas.ractive_;
  rclosed_ = cas.rclosed_;

  t2 = cas.t2all_[0]->at(0);
  den1  = cas.den1->clone();
  den2  = cas.den2->clone();

  rdm0_ = cas.rdm0all_->at(0);
  rdm1_ = cas.rdm1all_->at(0);
  rdm2_ = cas.rdm2all_->at(0);
  rdm3_ = cas.rdm3all_->at(0);
  rdm4_ = cas.rdm4all_->at(0);

  info_ = cas.info_;
  FCI_bare fci(info_->ciwfn());
  shared_ptr<const RDM<1>> ardm1;
  shared_ptr<const RDM<2>> ardm2;
  shared_ptr<const RDM<3>> ardm3;
  shared_ptr<const RDM<4>> ardm4;
  tie(ardm1, ardm2) = fci.rdm12_alpha(0,0);
  tie(ardm3, ardm4) = fci.rdm34_alpha(0,0);

  const int nclo = info_->nclosed();
  ardm1_ = fill_block<2,double>(ardm1, vector<int>(2,nclo), vector<IndexRange>(2,active_));
  ardm2_ = fill_block<4,double>(ardm2, vector<int>(4,nclo), vector<IndexRange>(4,active_));
  ardm3_ = fill_block<6,double>(ardm3, vector<int>(6,nclo), vector<IndexRange>(6,active_));
  ardm4_ = fill_block<8,double>(ardm4, vector<int>(8,nclo), vector<IndexRange>(8,active_));
}


void SPCASPT2::SPCASPT2::solve() {
  shared_ptr<Queue> dens2 = make_densityq();
  while (!dens2->done())
    dens2->next_compute();

  shared_ptr<Queue> dens1 = make_density1q();
  while (!dens1->done())
    dens1->next_compute();
}

#endif
