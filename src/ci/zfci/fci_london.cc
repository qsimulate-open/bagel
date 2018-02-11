//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_london.cc
// Copyright (C) 2017 Toru Shiozaki
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

#include <src/ci/zfci/fci_london.h>
#include <src/ci/zfci/jop_london.h>
#include <src/wfn/zreference.h>

using namespace std;
using namespace bagel;

FCI_London::FCI_London(shared_ptr<const PTree> a, shared_ptr<const Geometry> g, shared_ptr<const Reference> b,
                       const int ncore, const int nocc, shared_ptr<const ZCoeff_Block> coeff_zcas)
 : ZHarrison(a, g, b, ncore, nocc, coeff_zcas, false, false) {

  cout << "    * FCI using gauge including atomic orbitals" << endl;
  cout << "    * " << nele_ << " active electrons in " << norb_ << " orbitals."  << endl;

  auto rr = dynamic_pointer_cast<const ZReference>(ref_);
  if (!rr) throw runtime_error("FCI with GIAO requires an appropriate reference object");

  shared_ptr<const ZCoeff_Block> coeff = coeff_zcas;
  if (!coeff)
    coeff = init_coeff();
  update(coeff);
}


// obtain the coefficient matrix in striped format
shared_ptr<const ZCoeff_Block> FCI_London::init_coeff() {
  auto rr = dynamic_pointer_cast<const ZReference>(ref_);
  assert(rr);

  const int rndim = rr->zcoeff()->ndim();
  const int rmdim = rr->zcoeff()->mdim();
  ZMatrix coeff(rndim, 2*rmdim);
  for (int i = 0; i != rmdim; ++i) {
    copy_n(rr->zcoeff()->element_ptr(0,i), rndim, coeff.element_ptr(0,2*i));
    copy_n(rr->zcoeff()->element_ptr(0,i), rndim, coeff.element_ptr(0,2*i+1));
  }
  auto scoeff = make_shared<const ZCoeff_Striped>(coeff, ncore_, norb_, rr->nvirt());

  // Reorder as specified in the input so frontier orbitals contain the desired active space
  const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
  if (iactive) {
    set<int> active_indices;
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iactive)
      active_indices.insert(lexical_cast<int>(i->data()) - 1);
    scoeff = scoeff->set_active(active_indices, geom_->nele()-charge_);
  }
  return  scoeff->block_format();
}


void FCI_London::update(shared_ptr<const ZCoeff_Block> coeff) {
  Timer timer;
  jop_ = make_shared<Jop_London>(geom_, ncore_*2, (ncore_+norb_)*2, coeff);
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;
  const_denom();
}


