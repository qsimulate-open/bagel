//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relfci.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/ci/zfci/relfci.h>
#include <src/ci/zfci/reljop.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/mat1e/rel/reloverlap.h>
#include <src/util/exception.h>

using namespace std;
using namespace bagel;

RelFCI::RelFCI(shared_ptr<const PTree> a, shared_ptr<const Geometry> g, shared_ptr<const Reference> b,
               const int ncore, const int nocc, shared_ptr<const ZCoeff_Block> coeff_zcas, const bool store_c, const bool store_g)
 : ZHarrison(a, g, b, ncore, nocc, coeff_zcas, store_c, store_g) {

  cout << "    * Relativistic FCI" << endl;
  cout << "    * " << nele_ << " active electrons in " << norb_ << " orbitals."  << endl;
  cout << "    * gaunt    : " << (gaunt_ ? "true" : "false") << endl;
  cout << "    * breit    : " << (breit_ ? "true" : "false") << endl;

  auto rr = dynamic_pointer_cast<const RelReference>(ref_);
  if (!rr) throw runtime_error("RelFCI requires a relativistic reference object");

  gaunt_ = idata_->get<bool>("gaunt", rr->gaunt());
  breit_ = idata_->get<bool>("breit", rr->breit());
  assert(!geom_->magnetism());

  shared_ptr<const ZCoeff_Block> coeff = coeff_zcas;
  if (!coeff)
    coeff = init_coeff();
  update(coeff);

  // if integral dump is requested, do it here, and throw Termination
  const bool only_ints = idata_->get<bool>("only_ints", false);
  if (only_ints)
    dump_integrals_and_exit();
}


void RelFCI::dump_integrals_and_exit() const {
#ifndef DISABLE_SERIALIZATION
  auto rr = dynamic_pointer_cast<const RelReference>(ref_);
  assert(rr);
  OArchive ar("relref");
  auto rout = make_shared<RelReference>(geom_, jop_->coeff()->striped_format(), 0.0, rr->nneg(), rr->nclosed(), rr->nact(), rr->nvirt(), rr->gaunt(), rr->breit());
  ar << rout;
  dump_ints();
  throw Termination("Relativistic MO integrals are dumped on a file.");
#else
  throw runtime_error("You must compile with serialization in order to dump MO integrals into a file.");
#endif
}


// obtain the coefficient matrix in striped format
shared_ptr<const ZCoeff_Block> RelFCI::init_coeff() {
  // in case reference is not relativistic...
  if (!geom_->dfs())
    geom_ = geom_->relativistic(gaunt_);

  auto overlap = make_shared<RelOverlap>(geom_);
  auto hcore = make_shared<RelHcore>(geom_);

  auto rr = dynamic_pointer_cast<const RelReference>(ref_);
  assert(rr);

  shared_ptr<const ZCoeff_Striped> scoeff;
  if (rr->kramers()) {
    scoeff = rr->relcoeff_full();
  } else {
    // then compute Kramers adapated coefficient matrices
    scoeff = make_shared<const ZCoeff_Striped>(*rr->relcoeff_full(), ncore_, norb_, rr->relcoeff_full()->mdim()/4-ncore_-norb_, rr->relcoeff_full()->mdim()/2);
    scoeff = scoeff->init_kramers_coeff(geom_, overlap, hcore, 2*ref_->nclosed() + ref_->nact(), gaunt_, breit_);
  }

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


void RelFCI::update(shared_ptr<const ZCoeff_Block> coeff) {
  if (!geom_->dfs())
    geom_ = geom_->relativistic(gaunt_);

  Timer timer;
  jop_ = make_shared<RelJop>(geom_, ncore_*2, (ncore_+norb_)*2, coeff, gaunt_, breit_, store_half_ints_, store_gaunt_half_ints_);
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;
  const_denom();
}


