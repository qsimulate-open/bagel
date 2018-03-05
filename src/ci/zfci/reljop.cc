//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reljop.cc
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

#include <src/ci/zfci/reljop.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/mat1e/giao/relhcore_london.h>
#include <src/scf/dhf/dfock.h>

using namespace std;
using namespace bagel;


RelJop::RelJop(const shared_ptr<const Geometry> geom, const int nstart, const int nfence, shared_ptr<const ZCoeff_Block> coeff,
               const bool gaunt, const bool breit, const bool store_c, const bool store_g)
 : ZMOFile(geom, coeff), gaunt_(gaunt), breit_(breit) {
  init(nstart, nfence, store_c, store_g);
}


shared_ptr<ZMatrix> RelJop::compute_hcore() const {
  return geom_->magnetism() ? static_pointer_cast<ZMatrix>(make_shared<RelHcore_London>(geom_))
                            : static_pointer_cast<ZMatrix>(make_shared<RelHcore>(geom_));
}


shared_ptr<ZMatrix> RelJop::compute_fock(shared_ptr<const ZMatrix> hcore, const int nclosed, const bool store_c, const bool store_g) const {
  return make_shared<DFock>(geom_, hcore, coeff_->slice_copy(0, nclosed), gaunt_, breit_, store_c, /*robust*/breit_, /*scale_J*/1.0, /*scale_K*/1.0, store_g);
}


tuple<list<shared_ptr<RelDFHalf>>,list<shared_ptr<RelDFHalf>>>
  RelJop::compute_half(shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> coeff, const bool gaunt, const bool breit) {

  assert(!breit || gaunt);
  // (1) make DFDists
  vector<shared_ptr<const DFDist>> dfs;
  if (!gaunt) {
    dfs = geom->dfs()->split_blocks();
    dfs.push_back(geom->df());
  } else {
    dfs = geom->dfsl()->split_blocks();
  }
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, gaunt);

  // (2) first-transform
  list<shared_ptr<RelDFHalf>> half_complex = DFock::make_half_complex(dfdists, coeff);

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_exch, half_complex_exch2;
  for (auto& i : half_complex) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(/*docopy=*/false);
    i.reset();
    half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
    DFock::factorize(half_complex_exch);
  }
  half_complex.clear();

  if (breit) {
    // TODO Not the best implementation -- one could avoid apply_J to half-transformed objects
    auto breitint = make_shared<BreitInt>(geom);
    list<shared_ptr<Breit2Index>> breit_2index;
    for (int i = 0; i != breitint->Nblocks(); ++i) {
      breit_2index.push_back(make_shared<Breit2Index>(breitint->index(i), breitint->data(i), geom->df()->data2()));
      if (breitint->not_diagonal(i))
        breit_2index.push_back(breit_2index.back()->cross());
    }
    for (auto& i : half_complex_exch)
      half_complex_exch2.push_back(i->apply_J());

    for (auto& i : half_complex_exch)
      for (auto& j : breit_2index)
        if (i->alpha_matches(j)) {
          half_complex_exch2.push_back(i->apply_J()->multiply_breit2index(j));
          DFock::factorize(half_complex_exch2);
        }
  }
  return make_tuple(half_complex_exch, half_complex_exch2);
}


shared_ptr<Kramers<2,ZMatrix>> RelJop::compute_mo1e(shared_ptr<const Kramers<1,ZMatrix>> coeff) {
  auto out = make_shared<Kramers<2,ZMatrix>>();
  for (size_t i = 0; i != 4; ++i)
    out->emplace(i, make_shared<ZMatrix>(*coeff->at(i/2) % *core_fock_ * *coeff->at(i%2)));
  return out;
}


shared_ptr<Kramers<4,ZMatrix>> RelJop::compute_mo2e(shared_ptr<const Kramers<1,ZMatrix>> coeff) {

  auto compute = [&coeff, this](shared_ptr<Kramers<4,ZMatrix>> out, const bool gaunt, const bool breit) {
    array<list<shared_ptr<RelDFHalf>>,2> half_complex_exch, half_complex_exch2;
    for (int k = 0; k != 2; ++k)
      tie(half_complex_exch[k], half_complex_exch2[k]) = compute_half(geom_, coeff->at(k), gaunt, breit);

    // (4) compute (gamma|ii)
    auto full = make_shared<Kramers<2,ListRelDFFull>>();
    full->emplace({0,0}, compute_full(coeff->at(0), half_complex_exch[0], true));
    full->emplace({1,0}, compute_full(coeff->at(0), half_complex_exch[1], true));
    full->emplace({0,1}, compute_full(coeff->at(1), half_complex_exch[0], true));
    full->emplace({1,1}, compute_full(coeff->at(1), half_complex_exch[1], true));

    auto full2 = !breit ? full : make_shared<Kramers<2,ListRelDFFull>>();
    if (breit) {
      full2->emplace({0,0}, compute_full(coeff->at(0), half_complex_exch2[0], false));
      full2->emplace({1,0}, compute_full(coeff->at(0), half_complex_exch2[1], false));
      full2->emplace({0,1}, compute_full(coeff->at(1), half_complex_exch2[0], false));
      full2->emplace({1,1}, compute_full(coeff->at(1), half_complex_exch2[1], false));
    }

    // (5) compute 4-index quantities (16 of them - we are not using symmetry... and this is a very cheap step)
    const double gscale = gaunt ? (breit ? -0.25 : -1.0) : 1.0;
    for (size_t i = 0; i != 16; ++i) {
      // we do not need (1000, 0111, 1110, 0001, 1100, 0110)
      if (i == 8 || i == 7 || i == 14 || i == 1 || i == 12 || i == 6)
        continue;

      // we will construct (1010, 1101, 0100) later
      if (i == 10 || i == 13 || i == 4)
        continue;

      // we compute: 0000, 0010, 1001, 0101, 0011, 1011

      // TODO : put in if statement for apply_J if nact*nact is much smaller than number of MPI processes
      const int b2a = i/4;
      const int b2b = i%4;
      out->add(i, full->at(b2a)->form_4index(full2->at(b2b), gscale));

      // in breit cases we explicitly symmetrize the Hamiltnian (hence the prefactor 0.5 above)
      if (breit)
        *out->at(i) += *full2->at(b2a)->form_4index(full->at(b2b), gscale);
    }
  };

  auto out = make_shared<Kramers<4,ZMatrix>>();
  compute(out, false, false);

  if (gaunt_)
    compute(out, true, breit_);

  (*out)[{1,0,1,0}] = out->at({0,1,0,1})->clone();
  shared_ptr<ZMatrix> m1010 = out->at({0,1,0,1})->get_conjg();
  sort_indices<1,0,3,2,0,1,1,1>(m1010->data(), out->at({1,0,1,0})->data(), nocc_, nocc_, nocc_, nocc_);

  (*out)[{1,1,0,1}] = out->at({1,0,1,1})->clone();
  shared_ptr<ZMatrix> m1101 = out->at({1,0,1,1})->get_conjg();
  sort_indices<3,2,1,0,0,1,1,1>(m1101->data(), out->at({1,1,0,1})->data(), nocc_, nocc_, nocc_, nocc_);

  (*out)[{0,1,0,0}] = out->at({0,0,1,0})->clone();
  shared_ptr<ZMatrix> m0100 = out->at({0,0,1,0})->get_conjg();
  sort_indices<3,2,1,0,0,1,1,1>(m0100->data(), out->at({0,1,0,0})->data(), nocc_, nocc_, nocc_, nocc_);

  return out;
}


shared_ptr<ListRelDFFull> RelJop::compute_full(shared_ptr<const ZMatrix> coeff, list<shared_ptr<RelDFHalf>> half, const bool appj) {
  list<shared_ptr<const RelDFHalf>> halfc;
  for (auto& i : half)
    halfc.push_back(i);
  return compute_full(coeff, halfc, appj);
}


shared_ptr<ListRelDFFull> RelJop::compute_full(shared_ptr<const ZMatrix> coeff, list<shared_ptr<const RelDFHalf>> half, const bool appj) {
  // TODO remove once DFDistT class is fixed
  const bool transform_with_full = !(half.front()->nocc()*coeff->mdim() <= mpi__->size());
  if (!transform_with_full && appj) {
    for (auto& i : half)
      i = i->apply_J();
  }

  list<shared_ptr<RelDFFull>> dffull;
  for (auto& i : half)
    dffull.push_back(make_shared<RelDFFull>(i, coeff));
  DFock::factorize(dffull);

  assert(dffull.size() == 1 || dffull.size() == 3);

  for (auto& i : dffull)
    i->scale(i->fac());

  if (transform_with_full && appj)
    for (auto& i : dffull)
      i = i->apply_J();
  return make_shared<ListRelDFFull>(dffull);
}


