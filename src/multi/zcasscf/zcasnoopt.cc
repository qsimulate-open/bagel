//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcasnoopt.cc
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

#include <set>
#include <src/multi/zcasscf/zcasnoopt.h>
#include <src/prop/pseudospin/pseudospin.h>
#include <src/mat1e/rel/reloverlap.h>
#include <src/mat1e/rel/relhcore.h>

using namespace std;
using namespace bagel;


ZCASNoopt::ZCASNoopt(shared_ptr<const PTree> idat, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
 : ZCASNoopt_base(idat, geom, ref) {
  init();
  cout << "    * No orbital optimization will be performed!" << endl << endl;
}


ZCASNoopt_London::ZCASNoopt_London(shared_ptr<const PTree> idat, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
 : ZCASNoopt_base(idat, geom, ref) {
  init();
  cout << "    * No orbital optimization will be performed!" << endl << endl;
}


void ZCASNoopt::init_coeff() {
  auto relref = dynamic_pointer_cast<const RelReference>(ref_);
  auto scoeff = make_shared<const ZCoeff_Striped>(*relref->relcoeff_full(), nclosed_, nact_, nvirtnr_, nneg_);
  if (!relref->kramers()) {
    auto overlap = make_shared<RelOverlap>(geom_);
    auto hcore = make_shared<RelHcore>(geom_);
    scoeff = scoeff->init_kramers_coeff(geom_, overlap, hcore, 2*ref_->nclosed() + ref_->nact(), gaunt_, breit_);
  }
  coeff_ = scoeff->block_format();
}


void ZCASNoopt_London::init_coeff() {
  auto relref = dynamic_pointer_cast<const RelReference>(ref_);
  auto scoeff = make_shared<const ZCoeff_Striped>(*relref->relcoeff_full(), nclosed_, nact_, nvirtnr_, nneg_);
  coeff_ = scoeff->block_format();
}


void ZCASNoopt_base::compute() {
  if (external_rdm_.empty()) {
    cout << "    * Computing RDMs from FCI calculation " << endl;
    fci_->compute();
    fci_->compute_rdm12();
  } else {
    fci_->read_external_rdm12_av(external_rdm_);
  }
  energy_ = fci_->energy();

  if (canonical_)
    tie(coeff_,eig_,eigB_,occup_,occupB_) = semi_canonical_orb(kramers());

  // TODO When the Property class is implemented, this should be one
  shared_ptr<const PTree> aniso_data = idata_->get_child_optional("aniso");
  if (aniso_data) {
    if (geom_->magnetism()) {
      cout << "  ** Magnetic anisotropy analysis is currently only available for zero-field calculations; sorry." << endl;
    } else {
      const int nspin = aniso_data->get<int>("nspin", (idata_->get_vector<int>("state", 0)).size()-1);
      Pseudospin ps(nspin, geom_, fci_->conv_to_ciwfn(), aniso_data);
      ps.compute(energy_, coeff_->active_part());
    }
  }
}
