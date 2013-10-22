//
// BAGEL - Parallel electron correlation program.
// Filename: zharrison_rdm.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/zfci/zharrison.h>

using namespace std;
using namespace bagel;

void ZHarrison::compute_rdm12() {

  // first compute <00|00>
  for (auto& isp : space_->detmap()) { 
    shared_ptr<const Determinants> base_det = isp.second; 
    const int nelea = base_det->nelea();
    const int neleb = base_det->neleb();

    for (int istate = 0; istate != nstate_; ++istate) {
      shared_ptr<const ZCivec> cc = cc_->find(nelea, neleb)->data(istate);
      shared_ptr<const Determinants> int_det = int_space_->finddet(nelea-1, neleb-1);

      const int ij = norb_*norb_;
      auto d = make_shared<ZDvec>(int_det, ij);
      const int nri = d->size() / ij; 

      sigma_2e_annih_ab(cc, d);

      // (a^+b^+|ba), 
      auto rdm1 = make_shared<RDM<1,complex<double>>>(norb_);
      auto rdm2 = make_shared<RDM<2,complex<double>>>(norb_);
      assert(d->size() == ij*nri);

      zgemv_("c", nri, ij, 1.0, d->data(0)->data(), nri, cc->data(), 1, 0.0, rdm1->data(), 1);
      zgemm3m_("c", "N", ij, ij, nri, 1.0, d->data(0)->data(), nri, d->data(0)->data(), nri, 0.0, rdm2->data(), ij);
    }
  }

}
