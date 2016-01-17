//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: tardm.cc
// Copyright (C) 2015 Toru Shiozaki
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

#ifdef SPINFREEMETHOD_DETAIL

namespace {

template<int N>
using TA = TATensor<complex<double>,N>;
template<int N>
using MapType = map<pair<int, int>, shared_ptr<TA<N>>>;

template <int N>
inline vector<MapType<N+1>> annihilate_one(const vector<MapType<N>>& in, const IndexRange active, const int maxtile) {
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
        IndexRange ci;
        for (size_t i = 0; i != det->lena(); ++i)
          ci.merge(IndexRange("o", lenb, maxtile, ci.nblock(), ci.size()));

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

          // skip beta tiles
          if (lo[1] >= nact) continue;
          // target alpha string
          const bitset<nbit__> astring = det->string_bits_a(lo[0]/lenb);
          // number of loops and stride
          const size_t sstride = up[0]-lo[0];;
          const size_t tstride = (up[1]-lo[1]) * sstride;
          const size_t nloops = range.volume() / tstride;
          const size_t boff = lo[0] % lenb;
          for (int i = lo[1]; i != up[1]; ++i) {
            if (astring[i]) continue;
            const double sign = det->template sign<0>(astring, i);
            // source alpha string
            bitset<nbit__> s = astring; s.set(i);
            // obtain the ordinal of the source tile
            vector<size_t> lo2 {det_orig->template lexical<0>(s)*lenb + boff};
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
                  blas::ax_plus_y_n(sign, source.begin()+n*sstride, sstride, target.begin()+n*tstride+(i-lo[1])*sstride);
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
        IndexRange ci;
        for (size_t i = 0; i != det->lena(); ++i)
          ci.merge(IndexRange("o", lenb, maxtile, ci.nblock(), ci.size()));

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

          // skip alpha tiles
          if (lo[1] < nact) continue;

          // loop over TiledRange1 of the first index (column/row conversion, hence using back())
          for (auto& p : source->trange().data().back()) {
            if (p.first/lenb_orig != lo[0]/lenb)
              continue;
            vector<size_t> slo{p.first};
            vector<size_t> sup{p.second};
            for (int i = 2; i < lo.size(); ++i) {
              slo.push_back(lo[i]);
              sup.push_back(up[i]);
            }
            const vector<size_t> slo_r(slo.rbegin(), slo.rend());
            const size_t index = source->trange().tiles().ordinal(source->trange().element_to_tile(slo_r));

            const size_t sstride = sup[0]-slo[0];
            const size_t tstride = (up[1]-lo[1]) * (up[0]-lo[0]);
            const size_t t0stride = up[0]-lo[0];
            const size_t soff = slo[0] % lenb_orig;
            const size_t toff = lo[0] % lenb;
            const size_t volume = range.volume();

            taket->get_world().taskq.add(
              [=](typename TA<N+1>::value_type tile, typename TA<N>::value_type stile) {
                for (size_t ss = 0, tt = 0; tt != volume; ss += sstride, tt += tstride)
                  for (int i = lo[1]-nact; i != up[1]-nact; ++i)
                    for (auto& ts : det->string_bits_b()) {
                      if (ts[i]) continue;
                      bitset<nbit__> s = ts; s.set(i);
                      const size_t tlex = det->template lexical<1>(ts);
                      const size_t slex = det_orig->template lexical<1>(s);
                      if (tlex >= toff && tlex < toff+t0stride
                       && slex >= soff && slex < soff+sstride) {
                        const double sign = det->template sign<1>(ts, i);
                        tile[tt+tlex-toff+(i-lo[1]+nact)*t0stride] += sign*stile[ss+slex-soff];
                      }
                    }
              }, (*it).future(), source->find(index)
            );
          }
        }

        if (!exist)
          out[istate].emplace(cpair, taket);
      }
    }
  }
  return out;
}

} // end namespace


template<>
void SpinFreeMethod<complex<double>>::feed_rdm_ta() {

  // first set the CI vectors to TA tensors.
  shared_ptr<const RelCIWfn> ciwfn = info_->ciwfn();
  shared_ptr<const RelZDvec> reldvec = ciwfn->civectors();
  map<pair<int, int>, shared_ptr<const ZDvec>> dvecs = reldvec->dvecs();

  // nelea, neleb
  const int nele = dynamic_pointer_cast<const RelSpace>(reldvec->space())->nele();
  const int nstates = info_->ciwfn()->nstates();
  const int maxtile = info_->maxtile();

  vector<MapType<1>> tacivec(nstates);
  for (int istate = 0; istate != nstates; ++istate)
    for (auto& id : reldvec->dvecs()) {
      shared_ptr<const ZCivec> cc = id.second->data(istate);

      auto det = cc->det();
      IndexRange ci;
      for (size_t i = 0; i != det->lena(); ++i)
        ci.merge(IndexRange("o", det->lenb(), maxtile, ci.nblock(), ci.size()));
      auto taket = make_shared<TA<1>>({ci});

      for (auto it = taket->begin(); it != taket->end(); ++it) {
        auto range = taket->trange().make_tile_range(it.ordinal());
        auto lo = range.lobound();
        auto up = range.upbound();
        auto buf = make_shared<ZVectorB>(up[0]-lo[0]);
        copy_n(cc->data()+lo[0], up[0]-lo[0], buf->data());
        taket->init_tile(it, buf);
      }
      tacivec[istate].emplace(id.first, taket);
    }

  vector<MapType<2>> ta1vec = annihilate_one(tacivec, active_, maxtile);
  ta1vec[0].begin()->second->get_world().gop.fence();

  vector<MapType<3>> ta2vec;
  if (nele > 1) {
    ta2vec = annihilate_one(ta1vec, active_, maxtile);
    ta2vec[0].begin()->second->get_world().gop.fence();
  }

  vector<MapType<4>> ta3vec;
  if (nele > 2) {
    ta3vec = annihilate_one(ta2vec, active_, maxtile);
    ta3vec[0].begin()->second->get_world().gop.fence();
  }

  vector<MapType<5>> ta4vec;
  if (nele > 3) {
    ta4vec = annihilate_one(ta3vec, active_, maxtile);
    ta4vec[0].begin()->second->get_world().gop.fence();
  }

  // TODO rebuilding these quantities (remove the code in spinfreebase.cc once everything is done.
  rdm0all_ = make_shared<Vec<TATensor<complex<double>,0>>>();
  rdm1all_ = make_shared<Vec<TATensor<complex<double>,2>>>();
  rdm2all_ = make_shared<Vec<TATensor<complex<double>,4>>>();
  rdm3all_ = make_shared<Vec<TATensor<complex<double>,6>>>();
  rdm4all_ = make_shared<Vec<TATensor<complex<double>,8>>>();
  for (int ibra = 0; ibra != nstates; ++ibra)
    for (int iket = 0; iket != nstates; ++iket) { // TODO can be reduced by 2 using symmetry
      auto rdm0t = make_shared<TATensor<complex<double>,0>>(vector<IndexRange>());
      auto rdm1t = make_shared<TATensor<complex<double>,2>>(vector<IndexRange>(2,active_), true);
      auto rdm2t = make_shared<TATensor<complex<double>,4>>(vector<IndexRange>(4,active_), true);
      auto rdm3t = make_shared<TATensor<complex<double>,6>>(vector<IndexRange>(6,active_), true);
      auto rdm4t = make_shared<TATensor<complex<double>,8>>(vector<IndexRange>(8,active_), true);

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

      if (nele > 2)
        for (auto& i : ta3vec[ibra])
          for (auto& j : ta3vec[iket])
            if (i.first == j.first)
              (*rdm3t)("o4,o1,o5,o2,o6,o3") += (*i.second)("o0,o1,o2,o3").conj() * (*j.second)("o0,o4,o5,o6");

      if (nele > 3)
        for (auto& i : ta4vec[ibra])
          for (auto& j : ta4vec[iket])
            if (i.first == j.first)
              (*rdm4t)("o5,o1,o6,o2,o7,o3,o8,o4") += (*i.second)("o0,o1,o2,o3,o4").conj() * (*j.second)("o0,o5,o6,o7,o8");

      rdm0all_->emplace(ibra, iket, rdm0t);
      rdm1all_->emplace(ibra, iket, rdm1t);
      rdm2all_->emplace(ibra, iket, rdm2t);
      rdm3all_->emplace(ibra, iket, rdm3t);
      rdm4all_->emplace(ibra, iket, rdm4t);
    }
}

template<>
void SpinFreeMethod<double>::feed_rdm_ta() {
//assert(false);
}

#endif
