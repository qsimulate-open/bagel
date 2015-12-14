//
// BAGEL - Parallel electron correlation program.
// Filename: tardm.cc
// Copyright (C) 2015 Toru Shiozaki
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

#ifdef SPINFREEMETHOD_DETAIL

template<int N>
using TA = TATensor<complex<double>,N>;
template<int N>
using MapType = map<pair<int, int>, shared_ptr<TA<N>>>;


template<>
void SpinFreeMethod<complex<double>>::feed_rdm_ta() {

  // first set the CI vectors to TA tensors.
  shared_ptr<const RelCIWfn> ciwfn = info_->ciwfn();
  shared_ptr<const RelZDvec> reldvec = ciwfn->civectors();
  map<pair<int, int>, shared_ptr<const ZDvec>> dvecs = reldvec->dvecs();

  // nelea, neleb
  const int nele = dynamic_pointer_cast<const RelSpace>(reldvec->space())->nele();
  const int nstates = info_->ciwfn()->nstates();

  vector<MapType<1>> tacivec(nstates);

  for (int nket = 0; nket != nstates; ++nket) {
    for (auto& id : reldvec->dvecs()) {
      shared_ptr<const ZCivec> cket = id.second->data(nket);

      // index range is blocked by boundaries of beta strings
      IndexRange ci("o", cket->size(), cket->lenb());
      auto taket = make_shared<TA<1>>(vector<IndexRange>{ci});

      for (auto it = taket->begin(); it != taket->end(); ++it) {
        auto buf = make_shared<ZVectorB>(cket->lenb());
        auto range = taket->trange().make_tile_range(it.ordinal());
        auto lo = range.lobound();
        copy_n(cket->data()+lo[0], cket->lenb(), buf->data());
        taket->init_tile(it, buf);
      }

      tacivec[nket].emplace(id.first, taket);
    }
  }

  // checking norm
  assert(all_of(tacivec.begin(), tacivec.end(), [](MapType<1>& o) { return 1.0e-8 >
    fabs(1.0-accumulate(o.begin(), o.end(), 0.0, [](double n, pair<pair<int, int>, shared_ptr<TA<1>>> b) { return n+pow(b.second->norm(),2); })); }));

  // TODO take care of this later
  assert(nele > 1);

  auto space0 = dynamic_pointer_cast<const RelSpace>(reldvec->space());
  auto space1 = make_shared<RelSpace>(info_->nact(), nele-1);

  vector<MapType<2>> tadvec(nstates);
  for (int nket = 0; nket != nstates; ++nket) {
    for (auto& c : tacivec[nket]) {
      auto det_orig = space0->finddet(c.first.first, c.first.second);
      shared_ptr<const TA<1>> source = c.second;
      // annihilate a
      if (det_orig->nelea() > 0) {
        auto det = space1->finddet(det_orig->nelea()-1, det_orig->neleb());
        const size_t lenb = det->lenb();
        IndexRange ci("o", det->size(), lenb);

        const pair<int, int> cpair{det->nelea(), det->neleb()};
        const bool exist = tadvec[nket].find(cpair) != tadvec[nket].end();
        shared_ptr<TA<2>> taket = exist ? tadvec[nket].at(cpair) : make_shared<TA<2>>(vector<IndexRange>{ci, active_}, true);

        for (auto it = taket->begin(); it != taket->end(); ++it) {
          auto range = taket->trange().make_tile_range(it.ordinal());
          auto lo_r = range.lobound();
          auto up_r = range.upbound();
          const vector<size_t> lo(lo_r.rbegin(), lo_r.rend());
          const vector<size_t> up(up_r.rbegin(), up_r.rend());
          assert(up[0]-lo[0] == lenb && lo[0] % lenb == 0);

          // skip beta tiles
          if (lo[1] > info_->nact()) continue;
          // target alpha string
          const bitset<nbit__> astring = det->string_bits_a(lo[0]/lenb);
          for (int i = lo[1]; i != up[1]; ++i) {
            if (astring[i]) continue;
            const double sign = det->sign<0>(astring, i);
            // source alpha string
            bitset<nbit__> s = astring; s.set(i);
            // lexical order
            const size_t sa = det_orig->lexical<0>(s);
            // submit a task
            taket->get_world().taskq.add(
              [=](typename TA<2>::value_type target, typename TA<1>::value_type source) {
                blas::ax_plus_y_n(sign, source.begin(), lenb, target.begin()+(i-lo[1])*lenb);
              }, (*it).future(), source->find(sa)
            );
          }
        }
        if (!exist)
          tadvec[nket].emplace(cpair, taket);
      }

      // annihilate b
      if (det_orig->neleb() > 0) {
        auto det = space1->finddet(det_orig->nelea(), det_orig->neleb()-1);
        const size_t lenb = det->lenb();
        const size_t lenb_orig = det_orig->lenb();
        IndexRange ci("o", det->size(), det->lenb());

        const pair<int, int> cpair{det->nelea(), det->neleb()};
        const bool exist = tadvec[nket].find(cpair) != tadvec[nket].end();
        shared_ptr<TA<2>> taket = exist ? tadvec[nket].at(cpair) : make_shared<TA<2>>(vector<IndexRange>{ci, active_}, true);

        for (auto it = taket->begin(); it != taket->end(); ++it) {
          auto range = taket->trange().make_tile_range(it.ordinal());
          auto lo_r = range.lobound();
          auto up_r = range.upbound();
          const vector<size_t> lo(lo_r.rbegin(), lo_r.rend());
          const vector<size_t> up(up_r.rbegin(), up_r.rend());
          assert(up[0]-lo[0] == lenb && lo[0] % lenb == 0);

          // skip alpha tiles
          if (lo[1] <= info_->nact()) continue;

          const size_t sa = lo[0]/lenb;

          taket->get_world().taskq.add(
            [=](typename TA<2>::value_type target, typename TA<1>::value_type source) {
              for (int i = lo[1]; i != up[1]; ++i) {
                for (auto& ts : det->string_bits_b()) {
                  if (ts[i]) continue;
                  const double sign = det->sign<1>(ts, i);
                  bitset<nbit__> s = ts; s.set(i);
                  target[det->lexical<1>(s) + (i-lo[1])*lenb] += sign*source[det_orig->lexical<1>(s)];
                }
              }
            }, (*it).future(), source->find(sa)
          );
        }

        if (!exist)
          tadvec[nket].emplace(cpair, taket);
      }
    }
  }

  // TODO form 1RDM
#if 1
  set_rdm(0,0);
  shared_ptr<TA<2>> tmp = rdm1_->clone();
  tmp->fill_local(0.0);
  for (auto& i : tadvec[0])
    for (auto& j : tadvec[0])
      if (i.first == j.first)
        (*tmp)("o1,o2") += (*i.second)("o0,o1").conj() * (*i.second)("o0,o2");
  cout << *tmp << endl << "---" << endl << *rdm1_ << endl;
#endif
}

template<>
void SpinFreeMethod<double>::feed_rdm_ta() {
  assert(false);
}

#endif
