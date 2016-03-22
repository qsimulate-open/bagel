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

  // TODO compute and fill in ardms here
#if 0
  info_ = cas.info_;
  FCI_bare fci(info_->ciwfn());
#else
  ardm1_ = rdm1_;
  ardm2_ = rdm2_;
  ardm3_ = rdm3_;
  ardm4_ = rdm4_;
#endif
}


void SPCASPT2::SPCASPT2::solve() {
  shared_ptr<Queue> dens2 = make_densityq();
  while (!dens2->done())
    dens2->next_compute();

  den2->matrix()->print("second-order");

  shared_ptr<Queue> dens1 = make_density1q();
  while (!dens1->done())
    dens1->next_compute();

  den1->matrix()->print("first-order");
}


#endif
