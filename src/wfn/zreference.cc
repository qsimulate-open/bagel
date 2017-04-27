//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zreference.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/mat1e/mixedbasis.h>
#include <src/mat1e/giao/zoverlap.h>
#include <src/wfn/zreference.h>
#include <src/integral/compos/complexoverlapbatch.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::ZReference)

using namespace std;
using namespace bagel;

shared_ptr<Reference> ZReference::project_coeff(shared_ptr<const Geometry> geomin, const bool check_geom_change) const {
  assert(check_geom_change);
  if (!geomin->magnetism() || !geom_->magnetism())
    throw runtime_error("ZReference::project_coeff(...) is only for GIAO wavefunctions.  Not implemented for SOSCF or GIAO -> standard, etc.");

  shared_ptr<Reference> out;

  // standard 4-component wavefunction
  // project to a new basis
  const ZOverlap overlap(geomin);
  ZOverlap sinv = overlap;
  sinv.inverse();
  MixedBasis<ComplexOverlapBatch, ZMatrix> mixed(geom_, geomin);
  auto c = make_shared<ZCoeff>(sinv * mixed * *zcoeff_);

  // make coefficient orthogonal (under the overlap metric)
  ZMatrix unit = *c % overlap * *c;
  unit.inverse_half();
  *c *= unit;

  out = make_shared<ZReference>(geomin, c, energy_, nocc(), nact(), nvirt()+(geomin->nbasis()-geom_->nbasis()));

  return out;
}
