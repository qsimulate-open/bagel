//
// BAGEL - Parallel electron correlation program.
// Filename: ras/apply_operator.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/ras/apply_operator.h>
#include <src/math/algo.h>

using namespace std;
using namespace bagel;

namespace bagel {
struct HolesParticles {
  int ha, hb, pa, pb;
  HolesParticles() : ha(0), hb(0), pa(0), pb(0) {}
};
}

void ApplyOperator::operator()(const RASCivecView source, RASCivecView target, const vector<GammaSQ> operations, const vector<int> orbitals) const {
  assert(operations.size() == orbitals.size());

  const array<int, 3> ras = source.det()->ras();
  assert(ras == target.det()->ras());

  const int nops = operations.size();

  // First figure which block transformations are allowed, i.e., (h,p) -> (h',p') = (h+dh,p+dp)
  HolesParticles dhp;
  for (int op = 0; op < nops; ++op) {
    const int orb = orbitals[op];
    const bool spin = (operations[op]==GammaSQ::AnnihilateAlpha || operations[op]==GammaSQ::CreateAlpha);
    const int fac = (operations[op]==GammaSQ::CreateAlpha || operations[op]==GammaSQ::CreateBeta) ? 1 : -1;
    if (orb < ras[0]) // in RASI
      (spin ? dhp.ha : dhp.hb) -= fac;
    else if (orb >= ras[0]+ras[1]) // in RASIII
      (spin ? dhp.pa : dhp.pb) += fac;
    // else is in RASII and doesn't change hp
  }

  // loop over blocks in the target
  for (auto& tblock : target.blocks()) {
    if (!tblock) continue;
    // get appropriate block in source
    array<int, 4> tag = {tblock->stringsa()->nholes()     - dhp.ha, tblock->stringsb()->nholes()     - dhp.hb,
                       tblock->stringsa()->nparticles() - dhp.pa, tblock->stringsb()->nparticles() - dhp.pb};
    if (any_of(tag.begin(), tag.end(), [] (int i) { return i < 0; })) continue;

    auto sblock = source.block(tag[0], tag[1], tag[2], tag[3]);
    if (!sblock) continue;

    //const size_t sla = sblock->lena();
    const size_t slb = sblock->lenb();
    const size_t tla = tblock->lena();
    const size_t tlb = tblock->lenb();

    if (operations == vector<GammaSQ>{{GammaSQ::AnnihilateAlpha,GammaSQ::CreateAlpha}}) {
      shared_ptr<const RASDeterminants> det = target.det();
      assert(*det == *source.det());
      assert(tlb==slb);

      const int s = orbitals.front(); // creation operator on target
      const int r = orbitals.back();  // annihilation operator on target

      bitset<nbit__> mask1; mask1.set(r); mask1.set(s);
      bitset<nbit__> mask2; mask2.set(r);

      // (bit ^ maskrs) flips both r and s (if r!=s) or doesn't change (if r==s)
      bitset<nbit__> maskrs; maskrs.set(r); maskrs.flip(s);

      for (size_t ia = 0; ia < tla; ++ia) {
        auto tbit = tblock->string_bits_a(ia);
        if (((tbit & mask1) ^ mask2).none()) { // equivalent to tbit[s] && (r==s || !tbit[r])
          const bitset<nbit__> sbit = tbit ^ maskrs;
          const size_t slex = det->lexical_zero<0>(sbit);
          const int signrs = sign(sbit, r, s);
          blas::ax_plus_y_n(signrs, &sblock->element(0,slex), tlb, &tblock->element(0,ia));
        }
      }
    }
    else if (operations == vector<GammaSQ>{{GammaSQ::AnnihilateBeta,GammaSQ::CreateBeta}}) {
      shared_ptr<const RASDeterminants> det = target.det();
      assert(*det == *source.det());
      assert(tla==sblock->lena());

      const int s = orbitals.front(); // creation operator on target
      const int r = orbitals.back();  // annihilation operator on target

      bitset<nbit__> mask1; mask1.set(r); mask1.set(s);
      bitset<nbit__> mask2; mask2.set(r);

      // (bit ^ maskrs) flips both r and s (if r!=s) or doesn't change (if r==s)
      bitset<nbit__> maskrs; maskrs.set(r); maskrs.flip(s);

      for (size_t ib = 0; ib < tlb; ++ib) {
        auto tbit = tblock->string_bits_b(ib);
        if (((tbit & mask1) ^ mask2).none()) { // equivalent to tbit[s] && (r==s || !tbit[r])
          const bitset<nbit__> sbit = tbit ^ maskrs;
          const size_t slex = det->lexical_zero<1>(sbit);
          const int signrs = sign(sbit, r, s);
          daxpy_(tla, signrs, &sblock->element(slex,0), slb, &tblock->element(ib,0), tlb);
        }
      }
    }
    else
      throw logic_error("Not yet implemented!");
  }
}
