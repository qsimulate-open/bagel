//
// BAGEL - Parallel electron correlation program.
// Filename: relreference.cc
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

#include <src/rel/relreference.h>
#include <src/rel/reloverlap.h>
#include <src/integral/os/overlapbatch.h>
#include <src/integral/os/kineticbatch.h>

using namespace std;
using namespace bagel;

shared_ptr<Reference> RelReference::project_coeff(shared_ptr<const Geometry> geomin) const {
  // in this case we first form overlap matrices
  RelOverlap overlap(geomin);
  shared_ptr<ZMatrix> sinv = overlap.inverse();

  MixedBasis<OverlapBatch> smixed(geom_, geomin);
  MixedBasis<KineticBatch> tmixed(geom_, geomin);
  const int nb = geomin->nbasis();
  const int mb = geom_->nbasis();
  const complex<double> one(1.0);
  const complex<double> sca = one * (0.5/(c__*c__));
  ZMatrix mixed(nb*4, mb*4);
  mixed.copy_real_block(one,    0,    0, nb, mb, smixed);
  mixed.copy_real_block(one,   nb,   mb, nb, mb, smixed);
  mixed.copy_real_block(sca, 2*nb, 2*mb, nb, mb, tmixed);
  mixed.copy_real_block(sca, 3*nb, 3*mb, nb, mb, tmixed);

  auto c = make_shared<ZMatrix>(*sinv * mixed * *relcoeff_);
  return make_shared<RelReference>(geomin, c, energy_, 0, nocc(), nvirt()+2*(geomin->nbasis()-geom_->nbasis()), gaunt_, breit_);
}
