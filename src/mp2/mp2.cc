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
  const size_t naux = fullt->naux();

  cout << "    * 3-index integral transformation done" << endl;

  // make a list of static distribution of ij
  const int myrank = mpi__->rank();
  vector<vector<tuple<int,int,int,int>>> tasks(mpi__->size());
  {
    int nmax = 0;
    StaticDist ijdist(nocc*(nocc+1)/2, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0, cnt = 0; i < nocc; ++i)
        for (int j = i; j < nocc; ++j, ++cnt)
          if (cnt >= ijdist.start(inode) && cnt < ijdist.start(inode) + ijdist.size(inode))
            tasks[inode].push_back(make_tuple(j, i, /*mpitags*/-1,-1));
      if (tasks[inode].size() > nmax) nmax = tasks[inode].size();
    }
    for (auto& i : tasks) {
      const int n = i.size();
      for (int j = 0; j != nmax-n; ++j) i.push_back(make_tuple(-1,-1,-1,-1));
    }
  }
  const int nloop = tasks[0].size();

  // start communication (n fetch behind) - n is determined by memory size
  // the data is stored in a map
  map<int, shared_ptr<Matrix>> cache;
  // pair of node and set of integers
  vector<set<int>> cachetable(mpi__->size());

  auto cache_block = [&](const int nadd, const int ndrop) {
    assert(ndrop < nadd);
    if (ndrop >= 0) {
      for (int inode = 0; inode != mpi__->size(); ++inode) {
        const int id = get<0>(tasks[inode][ndrop]);
        const int jd = get<1>(tasks[inode][ndrop]);
        // if id and jd are no longer used in the cache, delete the element
        set<int> used;
        for (int i = ndrop+1; i <= nadd; ++i) {
          used.insert(get<0>(tasks[inode][i]));
          used.insert(get<1>(tasks[inode][i]));
        }
        if (!used.count(id)) {
          if (inode == myrank) cache.erase(id);
          cachetable[inode].erase(id);
        }
        if (!used.count(jd)) {
          if (inode == myrank) cache.erase(jd);
          cachetable[inode].erase(jd);
        }
      }
    }
    if (nadd < nloop) {
      // issue recv requests
      auto request_one_ = [&](const int i, const int rank) {
        if (i < 0) return -1;
        cachetable[rank].insert(i);
        int tag = -1;
        if (cache.find(i) == cache.end() && myrank == rank) {
          const int origin = fullt->locate(0, i*nvirt);
          if (origin == myrank) {
            cache[i] = fullt->get_slice(i*nvirt, (i+1)*nvirt).front();
          } else {
            cache[i] = make_shared<Matrix>(naux, nvirt, true);
            tag = mpi__->request_recv(cache[i]->data(), cache[i]->size(), origin, myrank*nocc+i);
          }
        }
        return tag;
      };

      // issue send requests
      auto send_one_ = [&](const int i, const int dest) {
        if (i < 0) return;
        // see if "i" is cached at dest
        if (cachetable[dest].count(i) || fullt->locate(0, i*nvirt) != myrank) return;
        mpi__->request_send(fullt->data() + (i*nvirt-fullt->bstart())*naux, nvirt*naux, dest, dest*nocc+i);
      };

      for (int inode = 0; inode != mpi__->size(); ++inode) {
        if (inode == myrank) {
          // recieve data from other processes
          get<2>(tasks[myrank][nadd]) = request_one_(get<0>(tasks[myrank][nadd]), myrank); // receive requests
          get<3>(tasks[myrank][nadd]) = request_one_(get<1>(tasks[myrank][nadd]), myrank);
        } else {
          // send data to other processes
          send_one_(get<0>(tasks[inode][nadd]), inode); // send requests
          send_one_(get<1>(tasks[inode][nadd]), inode);
          request_one_(get<0>(tasks[inode][nadd]), inode); // update cachetable
          request_one_(get<1>(tasks[inode][nadd]), inode);
        }
      }
    }
  };

  const int ncache = memory_size / (nvirt*nvirt*2);
  cout << "    * ncache = " << ncache << endl;
  for (int n = 0; n != min(ncache, nloop); ++n)
    cache_block(n, -1);

  // denominator info
  const vector<double> eig(ref_->eig().begin()+ncore_, ref_->eig().end());

  // loop over tasks
  energy_ = 0;
  for (int n = 0; n != nloop; ++n) {
    // take care of data. The communication should be hidden
    if (n+ncache < nloop)
      cache_block(n+ncache, n-1);

    const int i = get<0>(tasks[myrank][n]);
    const int j = get<1>(tasks[myrank][n]);
    if (i < 0 || j < 0) continue;

    const int ti = get<2>(tasks[myrank][n]);
    const int tj = get<3>(tasks[myrank][n]);
    if (ti >= 0) mpi__->wait(ti);
    if (tj >= 0) mpi__->wait(tj);

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


