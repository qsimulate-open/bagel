//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relfci_london.cc
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

#include <src/scf/dhf/dfock.h>
#include <src/ci/zfci/reljop.h>
#include <src/ci/zfci/relfci_london.h>
#include <src/mat1e/giao/relhcore_london.h>
#include <src/mat1e/giao/reloverlap_london.h>

using namespace std;
using namespace bagel;

RelFCI_London::RelFCI_London(shared_ptr<const PTree> a, shared_ptr<const Geometry> g, shared_ptr<const Reference> b,
                             const int ncore, const int nocc, shared_ptr<const ZCoeff_Block> coeff_zcas, const bool store_c, const bool store_g)
 : ZHarrison(a, g, b, ncore, nocc, coeff_zcas, store_c, store_g) {

  auto rr = dynamic_pointer_cast<const RelReference>(ref_);
  if (!rr) throw runtime_error("RelFCI_London requires a relativistic reference object");

  gaunt_ = idata_->get<bool>("gaunt", rr->gaunt());
  breit_ = idata_->get<bool>("breit", rr->breit());

  cout << "    * Relativistic FCI (with a magnetic field)" << endl;
  cout << "    * " << nele_ << " active electrons in " << norb_ << " orbitals."  << endl;
  cout << "    * gaunt    : " << (gaunt_ ? "true" : "false") << endl;
  cout << "    * breit    : " << (breit_ ? "true" : "false") << endl;

  // in case reference is not relativistic...
  if (!geom_->dfs())
    geom_ = geom_->relativistic(gaunt_);

  shared_ptr<const ZCoeff_Block> coeff = coeff_zcas;
  if (!coeff)
    coeff = init_coeff();
  update(coeff);

  // if integral dump is requested, do it here, and throw Termination
  const bool only_ints = idata_->get<bool>("only_ints", false);
  if (only_ints)
    dump_integrals_and_exit();
}


// obtain the coefficient matrix in striped format
shared_ptr<const ZCoeff_Block> RelFCI_London::init_coeff() {
  auto rr = dynamic_pointer_cast<const RelReference>(ref_);
  assert(rr);

  const int nvirt_nr = rr->relcoeff_full()->mdim()/4-ncore_-norb_;
  const int nneg = rr->relcoeff_full()->mdim()/2;
  auto scoeff = make_shared<const ZCoeff_Striped>(*rr->relcoeff_full(), ncore_, norb_, nvirt_nr, nneg, false);

  // Reorder as specified in the input so frontier orbitals contain the desired active space
  const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
  if (iactive) {
    set<int> active_indices;
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iactive)
      active_indices.insert(lexical_cast<int>(i->data()) - 1);
    scoeff = scoeff->set_active(active_indices, geom_->nele()-charge_);
  }
  return scoeff->block_format();
}


void RelFCI_London::update(shared_ptr<const ZCoeff_Block> coeff) {
  Timer timer;
  jop_ = make_shared<RelJop>(geom_, ncore_*2, (ncore_+norb_)*2, coeff, gaunt_, breit_, store_half_ints_, store_gaunt_half_ints_);
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;
  const_denom();
}


