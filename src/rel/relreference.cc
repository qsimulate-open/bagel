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
#include <src/molecule/zoverlap.h>
#include <src/london/reloverlap_london.h>
#include <src/integral/os/overlapbatch.h>
#include <src/integral/compos/complexoverlapbatch.h>
#include <src/integral/os/kineticbatch.h>
#include <src/integral/smallints1e_london.h>
#include <src/molecule/mixedbasis.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::RelReference)

using namespace std;
using namespace bagel;

shared_ptr<Reference> RelReference::project_coeff(shared_ptr<const Geometry> geomin) const {

  shared_ptr<Reference> out;
  const bool giao = (geomin->magnetism() || geom_->magnetism());

  // standard 4-component wavefunction
  if (rel_ && !giao) {
    // in this case we first form overlap matrices
    RelOverlap overlap(geomin);
    RelOverlap sinv = overlap;
    sinv.inverse();

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

    auto c = make_shared<ZMatrix>(sinv * mixed * *relcoeff_);

    // make coefficient orthogonal
    ZMatrix unit = *c % overlap * *c;
    unit.inverse_half();
    *c *= unit;

    out = make_shared<RelReference>(geomin, c, energy_, 0, nocc(), nvirt()+2*(geomin->nbasis()-geom_->nbasis()), gaunt_, breit_, rel_);

  // 4-component GIAO wavefunction
  } else if (rel_ && giao) {
    // in this case we first form overlap matrices
    RelOverlap_London overlap(geomin);
    RelOverlap_London sinv = overlap;
    sinv.inverse();

    shared_ptr<const Geometry> relgeomin = geomin->relativistic(false, false);
    MixedBasis<ComplexOverlapBatch, ZMatrix> smixed(geom_, geomin);
    MixedBasisArray<SmallInts1e_London<ComplexOverlapBatch>, ZMatrix> smallovl(geom_, relgeomin);

    const int nb = geomin->nbasis();
    const int mb = geom_->nbasis();
    const complex<double> r2 (0.25 / (c__*c__));
    const complex<double> i2 (0.0, r2.real());

    ZMatrix mixed(nb*4, mb*4);
    mixed.copy_block(0,    0, nb, mb, smixed);
    mixed.copy_block(nb,  mb, nb, mb, smixed);
    mixed.add_block( r2, 2*nb, 2*mb, nb, mb, *smallovl.data(0));
    mixed.add_block( r2, 3*nb, 3*mb, nb, mb, *smallovl.data(0));
    mixed.add_block( i2, 2*nb, 2*mb, nb, mb, *smallovl.data(1));
    mixed.add_block(-i2, 3*nb, 3*mb, nb, mb, *smallovl.data(1));
    mixed.add_block( i2, 2*nb, 3*mb, nb, mb, *smallovl.data(2));
    mixed.add_block( i2, 3*nb, 2*mb, nb, mb, *smallovl.data(2));
    mixed.add_block( r2, 2*nb, 3*mb, nb, mb, *smallovl.data(3));
    mixed.add_block(-r2, 3*nb, 2*mb, nb, mb, *smallovl.data(3));

    auto c = make_shared<ZMatrix>(sinv * mixed * *relcoeff_);

    // make coefficient orthogonal
    ZMatrix unit = *c % overlap * *c;
    unit.inverse_half();
    *c *= unit;

    out = make_shared<RelReference>(geomin, c, energy_, 0, nocc(), nvirt()+2*(geomin->nbasis()-geom_->nbasis()), gaunt_, breit_, rel_);

  // Non-relativistic GIAO wavefunction
  } else if (!rel_ && giao) {
    // project to a new basis
    const ZOverlap overlap(geomin);
    ZOverlap sinv = overlap;
    sinv.inverse();
    MixedBasis<ComplexOverlapBatch, ZMatrix> mixed(geom_, geomin);
    auto c = make_shared<ZCoeff>(sinv * mixed * *relcoeff_);

    // make coefficient orthogonal (under the overlap metric)
    ZMatrix unit = *c % overlap * *c;
    unit.inverse_half();
    *c *= unit;

    out = make_shared<RelReference>(geomin, c, energy_, 0, nocc(), nvirt()+(geomin->nbasis()-geom_->nbasis()), gaunt_, breit_, rel_);
    if (!geomin->magnetism())
      throw std::runtime_error("Projection from GIAO to real non-rel. basis would give complex coefficients.  Use the GIAO code at zero-field or restart.");

  } else {
    throw logic_error("Invalid RelReference formed");
  }
  return out;
}
