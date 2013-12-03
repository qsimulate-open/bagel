//
// BAGEL - Parallel electron correlation program.
// Filename: relreference.h
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

#ifndef __SRC_REL_RELREFERENCE_H
#define __SRC_REL_RELREFERENCE_H

#include <src/rel/relhcore.h>
#include <src/wfn/reference.h>
#include <src/molecule/mixedbasis.h>

namespace bagel {

class RelReference : public Reference {
  protected:
    const bool gaunt_;
    const bool breit_;
    const int nneg_;
    const std::shared_ptr<const ZMatrix> relcoeff_;
    const std::shared_ptr<const ZMatrix> relcoeff_full_;

  public:
    RelReference(std::shared_ptr<const Geometry> g, std::shared_ptr<const ZMatrix> c, const double en, const int nneg, const int nocc, const int nvirt, const bool ga, const bool br)
     : Reference(g, std::shared_ptr<Coeff>(), nocc, 0, nvirt, en), gaunt_(ga), breit_(br), nneg_(nneg), relcoeff_(c->slice(nneg_, c->mdim())), relcoeff_full_(c) {
    }

    const std::shared_ptr<const Coeff> coeff() const override { throw std::logic_error("RelReference::coeff() should not be called"); }
    const std::shared_ptr<const ZMatrix> relcoeff() const { return relcoeff_; }
    const std::shared_ptr<const ZMatrix> relcoeff_full() const { return relcoeff_full_; }

    bool gaunt() const { return gaunt_; }
    bool breit() const { return breit_; }
    int nneg() const { return nneg_; }

    std::shared_ptr<Reference> project_coeff(std::shared_ptr<const Geometry> geomin) const override;

};

}

#endif
