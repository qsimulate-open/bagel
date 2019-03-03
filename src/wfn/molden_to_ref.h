//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: molden_to_ref.h
// Copyright (C) 2018 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com> 
// Maintainer: NU theory
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

#ifndef __SRC_WFN_MOLDEN_TO_REF_H
#define __SRC_WFN_MOLDEN_TO_REF_H

#include <src/util/io/moldenin.h>
#include <src/wfn/zreference.h>
#include <src/wfn/relreference.h>

namespace bagel {
namespace {

std::shared_ptr<Reference> molden_to_ref(std::shared_ptr<const Geometry> geom, std::shared_ptr<const PTree> itree) {
  std::shared_ptr<Reference> ref;
  const std::string molden_file = itree->get<std::string>("molden_file");
  MoldenIn mfs(molden_file, geom->spherical());
  mfs.read();
  if (mfs.has_mo()) {
    MOInfo mo(geom);
    mfs >> mo;
    if (mo.coeff) {
      ref = std::make_shared<Reference>(geom);
      ref->set_coeff(mo.coeff);
      if (mo.coeffB)
        ref->set_coeff_AB(mo.coeff, mo.coeffB);
    } else if (mo.zcoeff) {
      auto zref = std::make_shared<ZReference>(geom);
      zref->set_zcoeff(mo.zcoeff);
      ref = zref;
    } else if (mo.relcoeff) {
      auto relref = std::make_shared<RelReference>(geom);
      relref->set_relcoeff(mo.relcoeff);
      ref = relref;
    } else {
      throw std::logic_error("MOInfo didn't have a MO. Strange");
    }
    ref->set_eig(mo.eig);
    ref->set_occup(mo.occup);
    if (mo.has_beta()) {
      ref->set_eigB(mo.eigB);
      ref->set_occupB(mo.occupB);
    }
  }
  return ref;
}

}
}

#endif
