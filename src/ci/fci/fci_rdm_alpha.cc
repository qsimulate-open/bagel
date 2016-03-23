//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_rdm_alpha.cc
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

#include <src/ci/fci/fci.h>

using namespace std;
using namespace bagel;

// This file is devoted for computing the so-called "alpha" density matrices
// <iaja>, <iaja,kl>, <iaja,kl,mn>, <iaja,kl,mn,op>
// the formulas are more or less the same.

tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> FCI::rdm12_alpha(const int ist, const int jst) {
  if (det_->compress()) {
    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
    cc_->set_det(detex);
  }

  shared_ptr<Civec> cbra = cc_->data(ist);
  shared_ptr<Civec> cket = cc_->data(jst);

  // iaja applied to bra
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  // kl applied to ket
  auto dket = make_shared<Dvec>(cket->det(), norb_*norb_);
  dket->zero();
  sigma_2a1(cket, dket);

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2, rdm2t;
  tie(rdm1,rdm2t) = compute_rdm12_last_step(dbra, dket, cbra);
  rdm2 = rdm2t->clone();
  blas::transpose(rdm2t->data(), norb_*norb_, norb_*norb_, rdm2->data());

  return tie(rdm1, rdm2);
}


tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>> FCI::rdm34_alpha(const int ist, const int jst) {
  shared_ptr<RDM<3>> rdm3;
  shared_ptr<RDM<4>> rdm4;
  return tie(rdm3, rdm4);
}
