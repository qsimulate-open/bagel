//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relreference.cc
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

#include <src/mat1e/mixedbasis.h>
#include <src/mat1e/rel/reloverlap.h>
#include <src/mat1e/giao/zoverlap.h>
#include <src/mat1e/giao/reloverlap_london.h>
#include <src/wfn/relreference.h>
#include <src/integral/os/overlapbatch.h>
#include <src/integral/compos/complexoverlapbatch.h>
#include <src/integral/os/kineticbatch.h>
#include <src/integral/smallints1e_london.h>
#include <src/ci/zfci/zharrison.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::RelReference)

using namespace std;
using namespace bagel;

shared_ptr<Reference> RelReference::project_coeff(shared_ptr<const Geometry> geomin, const bool check_geom_change) const {
  assert(check_geom_change);

  shared_ptr<Reference> out;
  const bool giao = (geomin->magnetism() || geom_->magnetism());

  // standard 4-component wavefunction
  if (!giao) {
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

    auto c2 = make_shared<RelCoeff_Striped>(*c, relcoeff_->nclosed(), relcoeff_->nact(), relcoeff_->nvirt_nr(), relcoeff_->nneg());
    out = make_shared<RelReference>(geomin, c2, energy_, nneg(), nocc(), nact(), nvirt()+2*(geomin->nbasis()-geom_->nbasis()), gaunt_, breit_, kramers_);

  // 4-component GIAO wavefunction
  } else {

    if (!geomin->magnetism() || !geom_->magnetism())
      throw std::runtime_error("Projection between GIAO and real basis sets is not implemented.   Use the GIAO code at zero-field or restart.");
    // in this case we first form overlap matrices
    RelOverlap_London overlap(geomin);
    RelOverlap_London sinv = overlap;
    sinv.inverse();

    shared_ptr<const Geometry> relgeomin = geomin->relativistic(false);
    MixedBasis<ComplexOverlapBatch, ZMatrix> smixed(geom_, relgeomin);
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

    auto c2 = make_shared<RelCoeff_Striped>(*c, relcoeff_->nclosed(), relcoeff_->nact(), relcoeff_->nvirt_nr(), relcoeff_->nneg());
    out = make_shared<RelReference>(geomin, c2, energy_, nneg(), nocc(), nact(), nvirt()+2*(geomin->nbasis()-geom_->nbasis()), gaunt_, breit_, kramers_);

  }
  return out;
}


shared_ptr<const Kramers<2,ZRDM<1>>> RelReference::rdm1(const int ist, const int jst) const {
  ZFCI_bare fci(ciwfn_);
  return fci.rdm1(ist, jst);
}


shared_ptr<const Kramers<4,ZRDM<2>>> RelReference::rdm2(const int ist, const int jst) const {
  ZFCI_bare fci(ciwfn_);
  return fci.rdm2(ist, jst);
}


shared_ptr<const Kramers<6,ZRDM<3>>> RelReference::rdm3(const int ist, const int jst) const {
  ZFCI_bare fci(ciwfn_);
  return fci.rdm3(ist, jst);
}


shared_ptr<const Kramers<8,ZRDM<4>>> RelReference::rdm4(const int ist, const int jst) const {
  ZFCI_bare fci(ciwfn_);
  return fci.rdm4(ist, jst);
}
