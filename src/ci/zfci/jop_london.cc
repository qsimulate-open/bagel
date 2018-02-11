//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: jop_london.cc
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

#include <src/ci/zfci/jop_london.h>
#include <src/mat1e/giao/zhcore.h>
#include <src/scf/giaohf/fock_london.h>

using namespace std;
using namespace bagel;


Jop_London::Jop_London(const shared_ptr<const Geometry> geom, const int nstart, const int nfence, shared_ptr<const ZCoeff_Block> coeff)
 : ZMOFile(geom, coeff) {
  init(nstart, nfence, false, false);
}


shared_ptr<ZMatrix> Jop_London::compute_hcore() const {
  return make_shared<ZHcore>(geom_);
}


shared_ptr<ZMatrix> Jop_London::compute_fock(shared_ptr<const ZMatrix> hcore, const int nclosed, const bool, const bool) const {
  auto ocoeff = coeff_->slice_copy(0, nclosed);
  ocoeff->scale(1.0/sqrt(2.0));
  return make_shared<Fock_London<1>>(geom_, hcore, nullptr, ocoeff, false, true);
}


shared_ptr<Kramers<2,ZMatrix>> Jop_London::compute_mo1e(shared_ptr<const Kramers<1,ZMatrix>> coeff) {
  auto out = make_shared<Kramers<2,ZMatrix>>();
  for (size_t i = 0; i != 4; ++i) {
    if (i/2 == i%2)
      out->emplace(i, make_shared<ZMatrix>(*coeff->at(i/2) % *core_fock_ * *coeff->at(i%2)));
    else
      out->emplace(i, out->at(0)->clone());
  }
  return out;
}


shared_ptr<Kramers<4,ZMatrix>> Jop_London::compute_mo2e(shared_ptr<const Kramers<1,ZMatrix>> coeff) {
  auto out = make_shared<Kramers<4,ZMatrix>>();

  auto df = dynamic_pointer_cast<const ComplexDFDist>(geom_->df());
  assert(df);

  array<shared_ptr<ComplexDFHalfDist>,2> half_complex_exch;
  for (int k = 0; k != 2; ++k)
    half_complex_exch[k] = df->complex_compute_half_transform(*coeff->at(k))->complex_apply_J();

  auto full = make_shared<Kramers<2,ComplexDFFullDist>>();
  full->emplace({0,0}, half_complex_exch[0]->complex_compute_second_transform(*coeff->at(0)));
  full->emplace({1,1}, half_complex_exch[1]->complex_compute_second_transform(*coeff->at(1)));

  // compute 4-index quantities
  for (size_t i = 0; i != 16; ++i) {
    // we do not need (1000, 0111, 1110, 0001, 1100, 0110)
    if (i == 8 || i == 7 || i == 14 || i == 1 || i == 12 || i == 6)
      continue;

    const int b2a = i/4;
    const int b2b = i%4;
    if (i == 0 || i == 3 || i == 12 || i == 15)
      out->add(i, full->at(b2a)->complex_form_4index(full->at(b2b), 1.0));
    else
      out->add(i, out->at(0)->clone());
  }

  return out;
}

