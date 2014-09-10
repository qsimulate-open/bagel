//
// BAGEL - Parallel electron correlation program.
// Filename: relreference_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/london/relreference_london.h>
#include <src/london/reloverlap_london.h>
#include <src/integral/compos/complexoverlapbatch.h>
#include <src/integral/smallints1e_london.h>
#include <src/molecule/mixedbasis.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::RelReference_London)

using namespace std;
using namespace bagel;

shared_ptr<Reference> RelReference_London::project_coeff(shared_ptr<const Geometry> geomin) const {
  // in this case we first form overlap matrices
  RelOverlap_London overlap(geomin);
  shared_ptr<ZMatrix> sinv = overlap.inverse();

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

  auto c = make_shared<ZMatrix>(*sinv * mixed * *relcoeff_);

  // make coefficient orthogonal
  ZMatrix unit = *c % overlap * *c;
  unit.inverse_half();
  *c *= unit;

  return make_shared<RelReference_London>(geomin, c, energy_, 0, nocc(), nvirt()+2*(geomin->nbasis()-geom_->nbasis()), gaunt_, breit_);
}
