//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zqvec.cc
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

#include <src/multi/zcasscf/zqvec.h>
#include <src/ci/zfci/reljop.h>

using namespace std;
using namespace bagel;


ZQvec::ZQvec(const int nbasis, const int nact, shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> rcoeff, shared_ptr<const ZMatrix> acoeff, const int nclosed,
             shared_ptr<const ZHarrison> fci, const bool gaunt, const bool breit)
 : ZMatrix(nbasis, nact*2) {

  assert(gaunt || !breit);
  assert(nbasis == rcoeff->mdim());

  auto compute = [&acoeff, &rcoeff, &geom, &fci, &nact] (const bool gaunt, const bool breit) {
    list<shared_ptr<RelDFHalf>> half, half2;
    tie(half, half2) = RelJop::compute_half(geom, acoeff, gaunt, breit);

    shared_ptr<ListRelDFFull> full = RelJop::compute_full(acoeff, half, true);
    shared_ptr<ListRelDFFull> full2 = !breit ? full : RelJop::compute_full(acoeff, half2, false);

    shared_ptr<ListRelDFFull> fullr = RelJop::compute_full(rcoeff, half, true);
    shared_ptr<ListRelDFFull> fullr2 = !breit ? fullr : RelJop::compute_full(rcoeff, half2, false);

    // form (rs|tu)*G(vs,tu) where r runs fastest
    shared_ptr<const ZMatrix> rdm2_av = fci->rdm2_av();

    auto full_d = make_shared<ListRelDFFull>();
    for (auto& ii : *full)
      full_d->push_back(ii->apply_2rdm(rdm2_av));

    shared_ptr<ListRelDFFull> full2_d = full_d;
    if (breit) {
      full2_d = make_shared<ListRelDFFull>();
      for (auto& ii : *full2)
        full2_d->push_back(ii->apply_2rdm(rdm2_av));
    }

    const double gscale = gaunt ? (breit ? -0.25 : -1.0) : 1.0;
    shared_ptr<ZMatrix> out = fullr->form_2index(full2_d, gscale, false);
    if (breit)
      *out += *fullr2->form_2index(full_d, gscale, false);
    return out;
  };

  *this = *compute(false, false);
  if (gaunt)
    *this += *compute(gaunt, breit);

  // complex conjugation due to bra-ket conjugation
  blas::conj_n(this->data(), this->size());
}
