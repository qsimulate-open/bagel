//
// BAGEL - Parallel electron correlation program.
// Filename: mp2.cc
// Copyright (C) 2012 Toru Shiozaki
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


// implements the MP2-F12 theory

#include <set>
#include <src/mp2/mp2.h>
#include <src/util/f77.h>
#include <src/scf/scf.h>
#include <src/smith/prim_op.h>
#include <src/mp2/f12int4.h>
#include <src/df/dfdistt.h>
#include <src/util/taskqueue.h>
#include <src/parallel/resources.h>

using namespace std;
using namespace bagel;

MP2::MP2(const shared_ptr<const PTree> input, const shared_ptr<const Geometry> g, const shared_ptr<const Reference> ref) : Method(input, g, ref) {

  scf_ = make_shared<SCF>(input, g, ref);
  scf_->compute();
  ref_ = scf_->conv_to_ref();

  cout << endl << "  === DF-MP2 calculation ===" << endl << endl;

  // checks for frozen core
  const bool frozen = idata_->get<bool>("frozen", false);
  ncore_ = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (ncore_) cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;

  if (geom_->df() == nullptr) throw logic_error("MP2 is only implemented with DF");

  // if three is a aux_basis keyword, we use that basis
  abasis_ = to_lower(idata_->get<string>("aux_basis", ""));

}


void MP2::compute() {
  const size_t nbasis = ref_->coeff()->mdim();
  const size_t nocc = ref_->nocc() - ncore_;
  if (nocc < 1) throw runtime_error("no correlated electrons");
  const size_t nvirt = nbasis - nocc - ncore_;
  if (nvirt < 1) throw runtime_error("no virtuals orbitals");


  shared_ptr<const Matrix> ocoeff = ref_->coeff()->slice(ncore_, ncore_+nocc);
  shared_ptr<const Matrix> vcoeff = ref_->coeff()->slice(ncore_+nocc, ncore_+nocc+nvirt);

  Timer timer;
  // compute transformed integrals
  shared_ptr<DFDistT> fullt;
  size_t memory_size;

  {
    // first compute half transformed integrals
    shared_ptr<DFHalfDist> half;
    if (abasis_.empty()) {
      half = geom_->df()->compute_half_transform(ocoeff);
    } else {
      auto info = make_shared<PTree>(); info->put("df_basis", abasis_);
      auto cgeom = make_shared<Geometry>(*geom_, info, false);
      half = cgeom->df()->compute_half_transform(ocoeff);
    }

    // second transform for virtual index and rearrange data
    {
      // this is now (naux, nvirt, nocc), distributed by nvirt*nocc. Always naux*nvirt block is localized to one node
      shared_ptr<DFFullDist> full = half->compute_second_transform(vcoeff)->apply_J()->swap();
      auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size(), full->nocc1());
      fullt = make_shared<DFDistT>(full, dist);
    }

    // the memory size info that was used for storing AO integrals here will be used for buffer
    memory_size = fullt->df()->block(0)->size();
    fullt->discard_df();
  }
  assert(fullt->nblocks() == 1);

  cout << "    * 3-index integral transformation done" << endl;

  // make a list of static distribution of ij
  const int myrank = mpi__->rank();
  vector<vector<tuple<int,int>>> tasks;

  {
    StaticDist ijdist(nocc*(nocc+1)/2, mpi__->size());
    int cnt = 0;
    for (int i = 0; i < nocc; ++i)
      for (int j = i; j < nocc; ++j, ++cnt)
        if (cnt >= ijdist.start(myrank) && cnt < ijdist.start(myrank) + ijdist.size(myrank))
          tasks.push_back(vector<tuple<int,int>>{ make_tuple(j, i) });
  }

  // start communication (n fetch behind) - n is determined by memory size
  // the data is stored in a map
  map<int, shared_ptr<const Matrix>> cache;

  auto cache_block = [&](const int nadd, const int ndrop) {
    assert(ndrop < nadd);
    if (nadd < tasks.size()) {
      const int ia = get<0>(tasks[nadd][myrank]);
      const int ja = get<1>(tasks[nadd][myrank]);
      if (cache.find(ia) == cache.end())
        cache[ia] = fullt->get_slice(ia*nvirt, (ia+1)*nvirt).front();
      if (cache.find(ja) == cache.end())
        cache[ja] = fullt->get_slice(ja*nvirt, (ja+1)*nvirt).front();
    }

    if (ndrop >= 0) {
      const int id = get<0>(tasks[ndrop][myrank]);
      const int jd = get<1>(tasks[ndrop][myrank]);
      // if id and jd are no longer used in the cache, delete the element
      set<int> used;
      for (int i = ndrop+1; i <= nadd; ++i) {
        used.insert(get<0>(tasks[i][myrank]));
        used.insert(get<1>(tasks[i][myrank]));
      }
      if (!used.count(id)) cache.erase(id);
      if (!used.count(jd)) cache.erase(jd);
    }
  };

cout << memory_size << endl;
  const int ncache = memory_size / (nvirt*nvirt*2);
  for (int n = 0; n != ncache; ++n)
    cache_block(n, -1);

  // denominator info
  const vector<double> eig(ref_->eig().begin()+ncore_, ref_->eig().end());

  // loop over tasks
  energy_ = 0;
  for (int n = 0; n != tasks.size(); ++n) {
    // take care of data
    const int m = n + ncache;
    const int dm = n - ncache;
    if (m < tasks.size())
      cache_block(m, dm);

    const int i = get<0>(tasks[n][myrank]);
    const int j = get<1>(tasks[n][myrank]);

    shared_ptr<const Matrix> iblock = cache.at(i);
    shared_ptr<const Matrix> jblock = cache.at(j);
    const Matrix mat(*iblock % *jblock);

    // should thread
    double en = 0.0;
    for (int a = 0; a != nvirt; ++a) {
      for (int b = a+1; b < nvirt; ++b) {
        const double ab = mat(a, b);
        const double ba = mat(b, a);
        en += 2.0*(ba*ba + ab*ab - ba*ab) / (-eig[a+nocc]+eig[i]-eig[b+nocc]+eig[j]);
      }
      const double aa = mat(a, a);
      en += aa*aa / (-eig[a+nocc]+eig[i]-eig[a+nocc]+eig[j]);
    }
    if (i != j) en *= 2.0;
    energy_ += en;
  }

  mpi__->allreduce(&energy_, 1);

  cout << "    * assembly done" << endl << endl;
  cout << "      MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << energy_ << setw(10) << setprecision(2) << timer.tick() << endl << endl;

  energy_ += ref_->energy();
  cout << "      MP2 total energy:       " << fixed << setw(15) << setprecision(10) << energy_ << endl << endl;

  // check if F12 is requested.
  const bool do_f12 = idata_->get<bool>("f12", false);
  if (do_f12) {
#ifdef HAVE_LIBSLATER
    const double gamma = idata_->get<double>("gamma", 1.5);
    cout << "    * F12 calculation requested with gamma = " << setprecision(2) << gamma << endl;
#if 0
    auto f12int = make_shared<F12Int>(idata_, geom_, ref_, gamma, ncore_);
#else
    auto f12ref = make_shared<F12Ref>(geom_, ref_, ncore_, gamma);
    f12ref->compute();
#endif
#else
  throw runtime_error("Slater-quadrature library not linked");
#endif
  }
}


