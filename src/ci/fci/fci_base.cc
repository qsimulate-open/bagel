//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_base.cc
// Copyright (C) 2016 Toru Shiozaki
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


#include <src/ci/fci/fci_base.h>

using namespace std;
using namespace bagel;


// note that this does not transform internal integrals (since it is not needed in CASSCF).
void FCI_base::rotate_rdms(shared_ptr<const Matrix> trans) {
  update_rdms(trans);
  jop_->update_1ext_ints(trans);
}


void FCI_base::update_rdms(shared_ptr<const Matrix> coeff) {
  for (auto& i : *rdm1_)
    i.second->transform(coeff);
  for (auto& i : *rdm2_)
    i.second->transform(coeff);

  // Only when #state > 1, this is needed.
  // Actually rdm1_av_ points to the same object as rdm1_ in 1 state runs. Therefore if you do twice, you get wrong.
  if (rdm1_->size() != 1) rdm1_av_->transform(coeff);
  if (rdm2_->size() != 1) rdm2_av_->transform(coeff);
  assert(rdm1_->size() != 1 || rdm1_->at(0) == rdm1_av_);
  assert(rdm2_->size() != 1 || rdm2_->at(0) == rdm2_av_);
}

