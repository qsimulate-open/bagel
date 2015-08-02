//
// BAGEL - Parallel electron correlation program.
// Filename: zqvec.cc
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

#include <src/multi/zcasscf/zqvec.h>
#include <src/scf/dhf/dfock.h>

using namespace std;
using namespace bagel;


ZQvec::ZQvec(const int nbasis, const int nact, shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> rcoeff, shared_ptr<const ZMatrix> acoeff, const int nclosed,
             shared_ptr<const ZHarrison> fci, const bool gaunt, const bool breit)
 : ZMatrix(nbasis, nact*2) {

  assert(gaunt || !breit);
  assert((*acoeff - *fci->jop()->coeff()).rms() < 1.0e-15);
  assert(nbasis == rcoeff->mdim());

  auto compute = [&acoeff, &rcoeff, &geom, &fci, &nact] (const bool gaunt, const bool breit) {
    list<shared_ptr<RelDFHalf>> half, half2;
    tie(half, half2) = RelMOFile::compute_half(geom, acoeff, gaunt, breit);

    list<shared_ptr<RelDFFull>> full = RelMOFile::compute_full(acoeff, half, true);
    list<shared_ptr<RelDFFull>> full2 = !breit ? full : RelMOFile::compute_full(acoeff, half2, false);

    list<shared_ptr<RelDFFull>> fullr = RelMOFile::compute_full(rcoeff, half, true);
    list<shared_ptr<RelDFFull>> fullr2 = !breit ? fullr : RelMOFile::compute_full(rcoeff, half2, false);

    // form (rs|tu)*G(vs,tu) where r runs fastest
    shared_ptr<const ZRDM<2>> rdm2_av = expand_kramers(fci->rdm2_av_kramers(), nact);

    list<shared_ptr<RelDFFull>> full_d;
    for (auto& ii : full)
      full_d.push_back(ii->apply_2rdm(rdm2_av));
    list<shared_ptr<RelDFFull>> full2_d;
    if (!breit)
      full2_d = full_d;
    else
      for (auto& ii : full2)
        full2_d.push_back(ii->apply_2rdm(rdm2_av));

    const double gscale = gaunt ? (breit ? -0.25 : -1.0) : 1.0;
    auto out = make_shared<ZMatrix>(rcoeff->mdim(), acoeff->mdim());
    {
      for (auto& ii : fullr)
        for (auto& jj : full2_d)
          if (ii->alpha_matches(jj))
            *out += *ii->form_2index(jj, gscale, false);
    }
    if (breit) {
      for (auto& ii : fullr2)
        for (auto& jj : full_d)
          if (ii->alpha_matches(jj))
            *out += *ii->form_2index(jj, gscale, false);
    }
    return out;
  };

  *this = *compute(false, false);
  if (gaunt)
    *this += *compute(gaunt, breit);
}
