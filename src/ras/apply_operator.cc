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

shared_ptr<const RASBlock<double>> ApplyOperator::get_block(const RASCivecView& source, array<int, 4>& dhp, shared_ptr<const RASBlock<double>> tblock) const {
  // get appropriate block in source
  array<int, 4> tag = {tblock->stringsa()->nholes()     - dhp[0], tblock->stringsb()->nholes()     - dhp[1],
                       tblock->stringsa()->nparticles() - dhp[2], tblock->stringsb()->nparticles() - dhp[3]};
  if (any_of(tag.begin(), tag.end(), [] (int i) { return i < 0; })) return nullptr;

  return source.block(tag[0], tag[1], tag[2], tag[3]);
}

void ApplyOperator::operator()(const double fac, const RASCivecView source, RASCivecView target, const vector<GammaSQ> operations, const vector<int> orbitals) const {
  assert(operations.size() == orbitals.size());

  shared_ptr<const RASDeterminants> sdet = source.det();
  shared_ptr<const RASDeterminants> tdet = target.det();

  const array<int, 3> ras = sdet->ras();
  assert(ras == tdet->ras());

  const int nops = operations.size();

  // First figure which block transformations are allowed, i.e., (h,p) -> (h',p') = (h+dh,p+dp)
  array<int, 4> dhp = {0, 0, 0, 0};
  for (int op = 0; op < nops; ++op) {
    const int orb = orbitals[op];
    const bool spin = (operations[op]==GammaSQ::AnnihilateAlpha || operations[op]==GammaSQ::CreateAlpha);
    const int delta = (operations[op]==GammaSQ::CreateAlpha || operations[op]==GammaSQ::CreateBeta) ? 1 : -1;
    if (orb < ras[0]) // in RASI
      (spin ? dhp[0] : dhp[1]) -= delta;
    else if (orb >= ras[0]+ras[1]) // in RASIII
      (spin ? dhp[2] : dhp[3]) += delta;
    // else is in RASII and doesn't change hp
  }

  if (operations == vector<GammaSQ>{{GammaSQ::AnnihilateAlpha,GammaSQ::CreateAlpha}}) {
    assert(*sdet == *tdet);

    const int r = orbitals.back();  // annihilation operator on target
    const int s = orbitals.front(); // creation operator on target

    bitset<nbit__> mask1; mask1.set(r); mask1.set(s);
    bitset<nbit__> mask2; mask2.set(r);

    // (bit ^ maskrs) flips both r and s (if r!=s) or doesn't change (if r==s)
    bitset<nbit__> maskrs; maskrs.set(r); maskrs.flip(s);

    // loop over blocks in the target
    for (auto& tblock : target.blocks()) {
      if (!tblock) continue;
      auto sblock = get_block(source, dhp, tblock);
      if (!sblock) continue;

      const size_t slb = sblock->lenb();
      const size_t tla = tblock->lena();
      const size_t tlb = tblock->lenb();

      for (size_t ia = 0; ia < tla; ++ia) {
        auto tbit = tblock->string_bits_a(ia);
        if (((tbit & mask1) ^ mask2).none()) { // equivalent to tbit[s] && (r==s || !tbit[r])
          const bitset<nbit__> sbit = tbit ^ maskrs;
          const size_t slex = sblock->stringsa()->lexical_zero(sbit);
          const int signrs = sign(sbit, r, s);
          blas::ax_plus_y_n(static_cast<double>(signrs)*fac, sblock->data() + slex*slb, tlb, tblock->data() + ia*tlb);
        }
      }
    }
  } else if (operations == vector<GammaSQ>{{GammaSQ::AnnihilateBeta,GammaSQ::CreateBeta}}) {
    assert(*sdet == *tdet);

    const int r = orbitals.back();  // annihilation operator on target
    const int s = orbitals.front(); // creation operator on target

    bitset<nbit__> mask1; mask1.set(r); mask1.set(s);
    bitset<nbit__> mask2; mask2.set(r);

    // (bit ^ maskrs) flips both r and s (if r!=s) or doesn't change (if r==s)
    bitset<nbit__> maskrs; maskrs.set(r); maskrs.flip(s);

    // loop over blocks in the target
    for (auto& tblock : target.blocks()) {
      if (!tblock) continue;
      auto sblock = get_block(source, dhp, tblock);
      if (!sblock) continue;

      const size_t slb = sblock->lenb();
      const size_t tla = tblock->lena();
      const size_t tlb = tblock->lenb();

      for (size_t ib = 0; ib < tlb; ++ib) {
        auto tbit = tblock->string_bits_b(ib);
        if (((tbit & mask1) ^ mask2).none()) { // equivalent to tbit[s] && (r==s || !tbit[r])
          const bitset<nbit__> sbit = tbit ^ maskrs;
          const size_t slex = sblock->stringsb()->lexical_zero(sbit);
          const int signrs = sign(sbit, r, s);
          daxpy_(tla, static_cast<double>(signrs)*fac, sblock->data() + slex, slb, tblock->data() + ib, tlb);
        }
      }
    }
  } else if (operations == vector<GammaSQ>{{GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta}}) {
    const int r = orbitals.back(); // beta annihilation on target
    const int s = orbitals.front();  // alpha creation on target

    // loop over blocks in the target
    for (auto& tblock : target.blocks()) {
      if (!tblock) continue;
      auto sblock = get_block(source, dhp, tblock);
      if (!sblock) continue;

      const size_t slb = sblock->lenb();
      const size_t tla = tblock->lena();
      const size_t tlb = tblock->lenb();

      // pre-process beta creation so it doesn't have to be done in a loop
      vector<tuple<size_t, int, size_t>> valid_bbits;
      const int nasign = (1 - ((tblock->stringsa()->nele()%2) << 1));
      for (size_t ib = 0; ib < tlb; ++ib) {
        bitset<nbit__> bbit = tblock->string_bits_b(ib);
        if (bbit[r])
          valid_bbits.emplace_back(ib, nasign*sign(bbit, r), sblock->stringsb()->lexical_zero(bbit ^ bitset<nbit__>(1 << r)));
      }

      if (!valid_bbits.empty()) {
        for (size_t ia = 0; ia < tla; ++ia) {
          bitset<nbit__> tabit = tblock->string_bits_a(ia);
          if (!tabit[s]) {
            bitset<nbit__> sabit = (tabit ^ bitset<nbit__>(1 << s));
            double* target_base = tblock->data() + ia*tlb;
            const double* source_base = sblock->data() + sblock->stringsa()->lexical_zero(sabit)*slb;
            const double c = fac * sign(sabit, s);
            for (auto& iex : valid_bbits)
              target_base[get<0>(iex)] += c * static_cast<double>(get<1>(iex)) * source_base[get<2>(iex)];
          }
        }
      }
    }
  } else if (operations == vector<GammaSQ>{{GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha}}) {
    const int r = orbitals.back(); // alpha annihilation on target
    const int s = orbitals.front();  // beta creation on target

    // loop over blocks in the target
    for (auto& tblock : target.blocks()) {
      if (!tblock) continue;
      auto sblock = get_block(source, dhp, tblock);
      if (!sblock) continue;

      const size_t slb = sblock->lenb();
      const size_t tla = tblock->lena();
      const size_t tlb = tblock->lenb();

      // pre-process beta creation so it doesn't have to be done in a loop
      vector<tuple<size_t, int, size_t>> valid_bbits;
      const int nasign = -(1 - ((tblock->stringsa()->nele()%2) << 1));
      for (size_t ib = 0; ib < tlb; ++ib) {
        bitset<nbit__> bbit = tblock->string_bits_b(ib);
        if (!bbit[s])
          valid_bbits.emplace_back(ib, nasign*sign(bbit, s), sblock->stringsb()->lexical_zero(bbit ^ bitset<nbit__>(1 << s)));
      }

      if (!valid_bbits.empty()) {
        for (size_t ia = 0; ia < tla; ++ia) {
          bitset<nbit__> tabit = tblock->string_bits_a(ia);
          if (tabit[r]) {
            bitset<nbit__> sabit = (tabit ^ bitset<nbit__>(1 << r));
            double* target_base = tblock->data() + ia*tlb;
            const double* source_base = sblock->data() + sblock->stringsa()->lexical_zero(sabit)*slb;
            const double c = fac * sign(sabit, r);
            for (auto& iex : valid_bbits)
              target_base[get<0>(iex)] += c * static_cast<double>(get<1>(iex)) * source_base[get<2>(iex)];
          }
        }
      }
    }
  }
  else
    throw logic_error("Not yet implemented!");
}
