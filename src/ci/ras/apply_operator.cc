//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/apply_operator.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/ci/ras/apply_block.h>
#include <src/ci/ras/apply_operator.h>
#include <src/util/math/algo.h>

using namespace std;
using namespace bagel;

shared_ptr<const RASBlock<double>> ApplyOperator::get_block(const RASCivecView& source, const array<int, 4>& dhp, const shared_ptr<const RASBlock<double>>& tblock) const {
  // get appropriate block in source
  array<int, 4> tag = {tblock->stringsa()->nholes()     - dhp[0], tblock->stringsb()->nholes()     - dhp[1],
                       tblock->stringsa()->nparticles() - dhp[2], tblock->stringsb()->nparticles() - dhp[3]};
  if (any_of(tag.begin(), tag.end(), [] (int i) { return i < 0; })) return nullptr;

  return source.block(tag[0], tag[1], tag[2], tag[3]);
}

shared_ptr<const CIBlockInfo<RASString>> ApplyOperator::get_blockinfo(const shared_ptr<const RASDeterminants>& sourcedet, const array<int, 4>& dhp, const shared_ptr<const CIBlockInfo<RASString>>& tblock) const {
  // get appropriate block in source
  array<int, 4> tag = {tblock->stringsa()->nholes()     - dhp[0], tblock->stringsb()->nholes()     - dhp[1],
                       tblock->stringsa()->nparticles() - dhp[2], tblock->stringsb()->nparticles() - dhp[3]};

  if (any_of(tag.begin(), tag.end(), [] (int i) { return i < 0; })) return nullptr;
  if (!sourcedet->allowed(tag[0], tag[1], tag[2], tag[3])) return nullptr;

  auto out = sourcedet->blockinfo(tag[0], tag[1], tag[2], tag[3]);
  return out->empty() ? nullptr : out;
}

void ApplyOperator::operator()(const double fac, const RASCivecView source, RASCivecView target, const vector<GammaSQ>& operations, const vector<int>& orbitals) const {
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
    const bool spin = is_alpha(operations[op]);
    const int delta = is_creation(operations[op]) ? 1 : -1;
    if (orb < ras[0]) // in RASI
      (spin ? dhp[0] : dhp[1]) -= delta;
    else if (orb >= ras[0]+ras[1]) // in RASIII
      (spin ? dhp[2] : dhp[3]) += delta;
    // else is in RASII and doesn't change hp
  }

  if (operations.size() == 1) {
    const int r = orbitals.front(); // operator
    const bool spin = is_alpha(operations.front());
    const bool action = is_creation(operations.front());
    RAS::Apply_block apply_block(r, action, spin);
    for (auto& tblock : target.blocks()) {
      if (!tblock) continue;
      auto sblock = get_block(source, dhp, tblock);
      if (sblock) apply_block(sblock, tblock, false);
    }
  } else if (operations == vector<GammaSQ>{{GammaSQ::CreateAlpha,GammaSQ::AnnihilateAlpha}}) {
    assert(*sdet == *tdet);

    const int r = orbitals.front();  // annihilation operator on target
    const int s = orbitals.back(); // creation operator on target

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
  } else if (operations == vector<GammaSQ>{{GammaSQ::CreateBeta,GammaSQ::AnnihilateBeta}}) {
    assert(*sdet == *tdet);

    const int r = orbitals.front();  // annihilation operator on target
    const int s = orbitals.back(); // creation operator on target

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
  } else if (operations == vector<GammaSQ>{{GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}}) {
    const int r = orbitals.front(); // beta annihilation on target
    const int s = orbitals.back();  // alpha creation on target

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
          valid_bbits.emplace_back(ib, nasign*sign(bbit, r), sblock->stringsb()->lexical_zero(bbit ^ (bitset<nbit__>(1) << r)));
      }

      if (!valid_bbits.empty()) {
        for (size_t ia = 0; ia < tla; ++ia) {
          bitset<nbit__> tabit = tblock->string_bits_a(ia);
          if (!tabit[s]) {
            bitset<nbit__> sabit = tabit ^ (bitset<nbit__>(1) << s);
            double* target_base = tblock->data() + ia*tlb;
            const double* source_base = sblock->data() + sblock->stringsa()->lexical_zero(sabit)*slb;
            const double c = fac * sign(sabit, s);
            for (auto& iex : valid_bbits)
              target_base[get<0>(iex)] += c * static_cast<double>(get<1>(iex)) * source_base[get<2>(iex)];
          }
        }
      }
    }
  } else if (operations == vector<GammaSQ>{{GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta}}) {
    const int r = orbitals.front(); // alpha annihilation on target
    const int s = orbitals.back();  // beta creation on target

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
          valid_bbits.emplace_back(ib, nasign*sign(bbit, s), sblock->stringsb()->lexical_zero(bbit ^ (bitset<nbit__>(1) << s)));
      }

      if (!valid_bbits.empty()) {
        for (size_t ia = 0; ia < tla; ++ia) {
          bitset<nbit__> tabit = tblock->string_bits_a(ia);
          if (tabit[r]) {
            bitset<nbit__> sabit = tabit ^ (bitset<nbit__>(1) << r);
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

void ApplyOperator::operator()(const double fac, const RASBlockVectors& source, RASBlockVectors& target, const vector<GammaSQ>& operations, const vector<int>& orbitals) const {
  assert(operations.size() == orbitals.size());

  const shared_ptr<const RASDeterminants>& sdet = source.det();
  const shared_ptr<const RASDeterminants>& tdet = target.det();

  const int M = source.mdim();
  if (M <= 0) return;

  const array<int, 3>& ras = sdet->ras();
  assert(ras == tdet->ras());

  const int nops = operations.size();

  // First figure which block transformations are allowed, i.e., (h,p) -> (h',p') = (h+dh,p+dp)
  array<int, 4> dhp = {0, 0, 0, 0};
  for (int op = 0; op < nops; ++op) {
    const int orb = orbitals[op];
    const bool spin = is_alpha(operations[op]);
    const int delta = is_creation(operations[op]) ? 1 : -1;
    if (orb < ras[0]) // in RASI
      (spin ? dhp[0] : dhp[1]) -= delta;
    else if (orb >= ras[0]+ras[1]) // in RASIII
      (spin ? dhp[2] : dhp[3]) += delta;
    // else is in RASII and doesn't change hp
  }

  const bool all_alpha = all_of(operations.begin(), operations.end(), [] (const GammaSQ& g) { return is_alpha(g); });
  const bool all_beta = all_of(operations.begin(), operations.end(), [] (const GammaSQ& g) { return is_beta(g); });

  assert(!(all_alpha && all_beta));

  if (all_alpha || all_beta) {
    vector<tuple<size_t, int, size_t, int, size_t, int>> base_daxpy;

    if (operations.size() == 1) {
      const int r = orbitals.front(); // operator
      const int base_phase = all_beta ? (1 - ((sdet->nelea()%2)<<1)) : 1;

      for (auto& tblock : tdet->blockinfo()) {
        if (tblock->empty()) continue;
        auto sblock = get_blockinfo(sdet, dhp, tblock);
        if (sblock) {
          const shared_ptr<const RASString>& source_strings = all_alpha ? sblock->stringsa() : sblock->stringsb();
          const shared_ptr<const RASString>& target_strings = all_alpha ? tblock->stringsa() : tblock->stringsb();

          const size_t tlen = target_strings->size();

          const size_t sstart = all_alpha ? sblock->lenb() : 1;
          const size_t tstart = all_alpha ? tblock->lenb() : 1;

          const size_t sstride = all_alpha ? 1 : sblock->lenb();
          const size_t tstride = all_alpha ? 1 : tblock->lenb();

          const size_t len = all_alpha ? sblock->lenb() : sblock->lena();
          assert(len == (all_alpha ? tblock->lenb() : tblock->lena()));

          for (size_t it = 0; it < tlen; ++it) {
            const bitset<nbit__> tbit = target_strings->strings(it);
            const bitset<nbit__> sbit = tbit ^ (bitset<nbit__>(1) << r);
            if (sbit.count() == source_strings->nele()) {
              assert(source_strings->contains(sbit));
              const size_t slex = source_strings->lexical_zero(sbit);
              const int rphase = sign(tbit, r);

              base_daxpy.emplace_back(len, base_phase*rphase, sblock->offset() + slex*sstart, sstride, tblock->offset() + it*tstart, tstride);
            }
          }
        }
      }
    } else if (operations == vector<GammaSQ>{{GammaSQ::CreateAlpha,GammaSQ::AnnihilateAlpha}}) {
      assert(*sdet == *tdet);

      const int r = orbitals.front();  // annihilation operator on target
      const int s = orbitals.back(); // creation operator on target

      bitset<nbit__> mask1; mask1.set(r); mask1.set(s);
      bitset<nbit__> mask2; mask2.set(r);

      // (bit ^ maskrs) flips both r and s (if r!=s) or doesn't change (if r==s)
      bitset<nbit__> maskrs; maskrs.set(r); maskrs.flip(s);

      // loop over blocks in the target
      for (auto& tblock : tdet->blockinfo()) {
        if (tblock->empty()) continue;
        auto sblock = get_blockinfo(sdet, dhp, tblock);
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
            base_daxpy.emplace_back(tlb, signrs, sblock->offset() + slex*slb, 1, tblock->offset() + ia*tlb, 1);
          }
        }
      }
    } else if (operations == vector<GammaSQ>{{GammaSQ::CreateBeta,GammaSQ::AnnihilateBeta}}) {
      assert(*sdet == *tdet);

      const int r = orbitals.front();  // annihilation operator on target
      const int s = orbitals.back(); // creation operator on target

      bitset<nbit__> mask1; mask1.set(r); mask1.set(s);
      bitset<nbit__> mask2; mask2.set(r);

      // (bit ^ maskrs) flips both r and s (if r!=s) or doesn't change (if r==s)
      bitset<nbit__> maskrs; maskrs.set(r); maskrs.flip(s);

      // loop over blocks in the target
      for (auto& tblock : tdet->blockinfo()) {
        if (tblock->empty()) continue;
        auto sblock = get_blockinfo(sdet, dhp, tblock);
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
            base_daxpy.emplace_back(tla, signrs, sblock->offset() + slex, slb, tblock->offset() + ib, tlb);
          }
        }
      }
    }
    else {
      throw logic_error("Not yet implemented!");
    }

    for (int m = 0; m < M; ++m) {
      for (auto& i : base_daxpy) {
        const double f = static_cast<double>(get<1>(i)) * fac;
        const double* source_data = source.element_ptr(get<2>(i), m);
        double* target_data = target.element_ptr(get<4>(i), m);
        daxpy_(get<0>(i), f, source_data, get<3>(i), target_data, get<5>(i));
      }
    }
  } else {
    unordered_map<size_t, vector<tuple<size_t, int, size_t>>> abits_map;
    unordered_map<size_t, vector<tuple<size_t, int, size_t>>> bbits_map;

    // <stringatag, stringbtag, source offset, source lenb, target offset, target lenb
    vector<array<size_t, 6>> block_info;

    if (operations == vector<GammaSQ>{{GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}}) {
      const int r = orbitals.front(); // beta annihilation on target
      const int s = orbitals.back();  // alpha creation on target

      // loop over blocks in the target
      for (auto& tblock : tdet->blockinfo()) {
        if (tblock->empty()) continue;
        auto sblock = get_blockinfo(sdet, dhp, tblock);
        if (sblock) {
          const size_t slb = sblock->lenb();
          const size_t tla = tblock->lena();
          const size_t tlb = tblock->lenb();

          bool compute_block = true;

          if (bbits_map.find(tblock->stringsb()->tag())==bbits_map.end()) {
            vector<tuple<size_t, int, size_t>> valid_bbits;
            const int nasign = (1 - ((tblock->stringsa()->nele()%2) << 1));
            for (size_t ib = 0; ib < tlb; ++ib) {
              const bitset<nbit__> bbit = tblock->string_bits_b(ib);
              if (bbit[r])
                valid_bbits.emplace_back(ib, nasign*sign(bbit, r), sblock->stringsb()->lexical_zero(bbit ^ (bitset<nbit__>(1) << r)));
            }
            if (!valid_bbits.empty())
              bbits_map.emplace(tblock->stringsb()->tag(), move(valid_bbits));
            else
              compute_block = false;
          }

          if (abits_map.find(tblock->stringsa()->tag())==abits_map.end()) {
            vector<tuple<size_t, int, size_t>> valid_abits;
            for (size_t ia = 0; ia < tla; ++ia) {
              const bitset<nbit__> abit = tblock->string_bits_a(ia);
              if (!abit[s])
                valid_abits.emplace_back(ia, sign(abit, s), sblock->stringsa()->lexical_zero(abit ^ (bitset<nbit__>(1) << s)));
            }
            if (!valid_abits.empty())
              abits_map.emplace(tblock->stringsa()->tag(), move(valid_abits));
            else
              compute_block = false;
          }

          if (compute_block)
            block_info.emplace_back(array<size_t, 6>{{tblock->stringsa()->tag(), tblock->stringsb()->tag(),
                                     sblock->offset(), slb, tblock->offset(), tlb}});
        }
      }
    } else if (operations == vector<GammaSQ>{{GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta}}) {
      const int r = orbitals.front(); // alpha annihilation on target
      const int s = orbitals.back();  // beta creation on target

      // loop over blocks in the target
      for (auto& tblock : tdet->blockinfo()) {
        if (tblock->empty()) continue;
        auto sblock = get_blockinfo(sdet, dhp, tblock);
        if (sblock) {
          const size_t slb = sblock->lenb();
          const size_t tla = tblock->lena();
          const size_t tlb = tblock->lenb();

          bool compute_block = true;

          if (bbits_map.find(tblock->stringsb()->tag())==bbits_map.end()) {
            vector<tuple<size_t, int, size_t>> valid_bbits;
            const int nasign = -(1 - ((tblock->stringsa()->nele()%2) << 1));
            for (size_t ib = 0; ib < tlb; ++ib) {
              bitset<nbit__> bbit = tblock->string_bits_b(ib);
              if (!bbit[s])
                valid_bbits.emplace_back(ib, nasign*sign(bbit, s), sblock->stringsb()->lexical_zero(bbit ^ (bitset<nbit__>(1) << s)));
            }
            if (!valid_bbits.empty())
              bbits_map.emplace(tblock->stringsb()->tag(), move(valid_bbits));
            else
              compute_block = false;
          }

          if (abits_map.find(tblock->stringsa()->tag())==abits_map.end()) {
            vector<tuple<size_t, int, size_t>> valid_abits;
            for (size_t ia = 0; ia < tla; ++ia) {
              const bitset<nbit__> abit = tblock->string_bits_a(ia);
              if (abit[r])
                valid_abits.emplace_back(ia, sign(abit, r), sblock->stringsa()->lexical_zero(abit ^ (bitset<nbit__>(1) << r)));
            }
            if (!valid_abits.empty())
              abits_map.emplace(tblock->stringsa()->tag(), move(valid_abits));
            else
              compute_block = false;
          }

          if (compute_block)
            block_info.emplace_back(array<size_t, 6>{{tblock->stringsa()->tag(), tblock->stringsb()->tag(),
                                     sblock->offset(), slb, tblock->offset(), tlb}});
        }
      }
    }
    else
      throw logic_error("Not yet implemented!");

    for (int m = 0; m < M; ++m) {
      for (auto& block : block_info) {
        for (auto& a : abits_map[block[0]]) {
          const double* source_base = source.element_ptr(block[2] + block[3]*get<2>(a), m);
          double* target_base = target.element_ptr(block[4] + block[5]*get<0>(a), m);
          const double f = fac * static_cast<double>(get<1>(a));
          for (auto& b : bbits_map[block[1]]) {
            target_base[get<0>(b)] += f * static_cast<double>(get<1>(b)) * source_base[get<2>(b)];
          }
        }
      }
    }
  }
}
