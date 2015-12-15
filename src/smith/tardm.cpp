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

template <int N>
inline vector<MapType<N+1>> annihilate_one(const vector<MapType<N>>& in, const IndexRange active) {
  const int nele = in[0].begin()->first.first + in[0].begin()->first.second;
  const int nact = active.size()/2;
  const int nstates = in.size();
  assert(nact*2 == active.size());
  auto ssp = make_shared<RelSpace>(nact, nele);
  auto tsp = make_shared<RelSpace>(nact, nele-1);

  vector<MapType<N+1>> out(nstates);
  for (int istate = 0; istate != nstates; ++istate) {
    for (auto& c : in[istate]) {
      auto det_orig = ssp->finddet(c.first.first, c.first.second);
      shared_ptr<const TA<N>> source = c.second;
      // annihilate a
      if (det_orig->nelea() > 0) {
        auto det = tsp->finddet(det_orig->nelea()-1, det_orig->neleb());
        const size_t lenb = det->lenb();
        IndexRange ci("o", det->size(), lenb);

        const pair<int, int> cpair{det->nelea(), det->neleb()};
        const bool exist = out[istate].find(cpair) != out[istate].end();
        vector<IndexRange> ind{ci};
        for (int i = 0; i != N; ++i)
          ind.push_back(active);
        shared_ptr<TA<N+1>> taket = exist ? out[istate].at(cpair) : make_shared<TA<N+1>>(ind, true);

        for (auto it = taket->begin(); it != taket->end(); ++it) {
          auto range = taket->trange().make_tile_range(it.ordinal());
          auto lo_r = range.lobound();
          auto up_r = range.upbound();
          const vector<size_t> lo(lo_r.rbegin(), lo_r.rend());
          const vector<size_t> up(up_r.rbegin(), up_r.rend());
          assert(up[0]-lo[0] == lenb && lo[0] % lenb == 0);

          // skip beta tiles
          if (lo[1] >= nact) continue;
          // target alpha string
          const bitset<nbit__> astring = det->string_bits_a(lo[0]/lenb);
          // number of loops and stride
          const size_t sstride = lenb;
          const size_t tstride = (up[1]-lo[1]) * lenb;
          const size_t nloops = range.volume() / tstride;
          for (int i = lo[1]; i != up[1]; ++i) {
            if (astring[i]) continue;
            const double sign = det->template sign<0>(astring, i);
            // source alpha string
            bitset<nbit__> s = astring; s.set(i);
            // obtain the ordinal of the source tile
            vector<size_t> lo2(1, det_orig->template lexical<0>(s)*lenb);
            for (int i = 2; i < lo.size(); ++i)
              lo2.push_back(lo[i]);
            const vector<size_t> lo2_r(lo2.rbegin(), lo2.rend());
            const size_t sa = source->trange().tiles().ordinal(source->trange().element_to_tile(lo2_r));

            // submit a task
            taket->get_world().taskq.add(
              [=](typename TA<N+1>::value_type target, typename TA<N>::value_type source) {
                assert(target.size() == tstride*nloops);
                assert(source.size() == sstride*nloops);
                for (size_t n = 0; n != nloops; ++n)
                  blas::ax_plus_y_n(sign, source.begin()+n*sstride, lenb, target.begin()+n*tstride+(i-lo[1])*lenb);
              }, (*it).future(), source->find(sa)
            );
          }
        }
        if (!exist)
          out[istate].emplace(cpair, taket);
      }

      // annihilate b
      if (det_orig->neleb() > 0) {
        auto det = tsp->finddet(det_orig->nelea(), det_orig->neleb()-1);
        const size_t lenb = det->lenb();
        const size_t lenb_orig = det_orig->lenb();
        IndexRange ci("o", det->size(), lenb);

        const pair<int, int> cpair{det->nelea(), det->neleb()};
        const bool exist = out[istate].find(cpair) != out[istate].end();
        vector<IndexRange> ind{ci};
        for (int i = 0; i != N; ++i)
          ind.push_back(active);
        shared_ptr<TA<N+1>> taket = exist ? out[istate].at(cpair) : make_shared<TA<N+1>>(ind, true);

        for (auto it = taket->begin(); it != taket->end(); ++it) {
          auto range = taket->trange().make_tile_range(it.ordinal());
          auto lo_r = range.lobound();
          auto up_r = range.upbound();
          const vector<size_t> lo(lo_r.rbegin(), lo_r.rend());
          const vector<size_t> up(up_r.rbegin(), up_r.rend());
          assert(up[0]-lo[0] == lenb && lo[0] % lenb == 0);

          // skip alpha tiles
          if (lo[1] < nact) continue;

          vector<size_t> lo2(1, lo[0]/lenb*lenb_orig);
          for (int i = 2; i < lo.size(); ++i)
            lo2.push_back(lo[i]);
          const vector<size_t> lo2_r(lo2.rbegin(), lo2.rend());
          const size_t sa = source->trange().tiles().ordinal(source->trange().element_to_tile(lo2_r));

          const size_t sstride = lenb_orig;
          const size_t tstride = (up[1]-lo[1]) * lenb;
          const size_t nloops = range.volume() / tstride;

          taket->get_world().taskq.add(
            [=](typename TA<N+1>::value_type tile, typename TA<N>::value_type stile) {
              for (size_t n = 0, soff = 0, toff = 0; n != nloops; ++n, soff += sstride, toff += tstride)
                for (int i = lo[1]-nact; i != up[1]-nact; ++i)
                  for (auto& ts : det->string_bits_b()) {
                    if (ts[i]) continue;
                    const double sign = det->template sign<1>(ts, i);
                    bitset<nbit__> s = ts; s.set(i);
                    tile[toff+det->template lexical<1>(ts)+(i-lo[1]+nact)*lenb] += sign*stile[soff+det_orig->template lexical<1>(s)];
                  }
            }, (*it).future(), source->find(sa)
          );
        }

        if (!exist)
          out[istate].emplace(cpair, taket);
      }
    }
  }
  return out;
}


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
  for (int istate = 0; istate != nstates; ++istate)
    for (auto& id : reldvec->dvecs()) {
      shared_ptr<const ZCivec> cc = id.second->data(istate);

      // index range is blocked by boundaries of beta strings
      IndexRange ci("o", cc->size(), cc->lenb());
      auto taket = make_shared<TA<1>>(vector<IndexRange>{ci});

      for (auto it = taket->begin(); it != taket->end(); ++it) {
        auto buf = make_shared<ZVectorB>(cc->lenb());
        auto range = taket->trange().make_tile_range(it.ordinal());
        auto lo = range.lobound();
        copy_n(cc->data()+lo[0], cc->lenb(), buf->data());
        taket->init_tile(it, buf);
      }
      tacivec[istate].emplace(id.first, taket);
    }

  vector<MapType<2>> ta1vec = annihilate_one(tacivec, active_);

// TODO understand the reason why we need this!! //////
ta1vec[0].begin()->second->get_world().gop.fence();
///////////////////////////////////////////////////////

  vector<MapType<3>> ta2vec;
  if (nele > 1)
    ta2vec = annihilate_one(ta1vec, active_);

  // TODO rebuilding these quantities (remove the code in spinfreebase.cc once everything is done.
  rdm0all_ = make_shared<Vec<TATensor<complex<double>,0>>>();
  rdm1all_ = make_shared<Vec<TATensor<complex<double>,2>>>();
  rdm2all_ = make_shared<Vec<TATensor<complex<double>,4>>>();
  for (int ibra = 0; ibra != nstates; ++ibra)
    for (int iket = 0; iket != nstates; ++iket) { // TODO can be reduced by 2 using symmetry
      auto rdm0t = make_shared<TATensor<complex<double>,0>>(vector<IndexRange>());
      auto rdm1t = make_shared<TATensor<complex<double>,2>>(vector<IndexRange>(2,active_), true);
      auto rdm2t = make_shared<TATensor<complex<double>,4>>(vector<IndexRange>(4,active_), true);

      (*rdm0t)("") = ibra == iket ? 1.0 : 0.0;

      for (auto& i : ta1vec[ibra])
        for (auto& j : ta1vec[iket])
          if (i.first == j.first)
            (*rdm1t)("o2,o1") += (*i.second)("o0,o1").conj() * (*j.second)("o0,o2"); // note: reversing indices

      if (nele > 1)
        for (auto& i : ta2vec[ibra])
          for (auto& j : ta2vec[iket])
            if (i.first == j.first)
              (*rdm2t)("o3,o1,o4,o2") += (*i.second)("o0,o1,o2").conj() * (*j.second)("o0,o3,o4");

      rdm0all_->emplace(ibra, iket, rdm0t);
      rdm1all_->emplace(ibra, iket, rdm1t);
      rdm2all_->emplace(ibra, iket, rdm2t);
    }
}

template<>
void SpinFreeMethod<double>::feed_rdm_ta() {
  assert(false);
}

#endif
