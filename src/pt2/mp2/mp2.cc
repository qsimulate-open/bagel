//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mp2.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/scf/hf/rhf.h>
#include <src/df/dfdistt.h>
#include <src/pt2/mp2/mp2.h>
#include <src/pt2/mp2/mp2cache.h>
#include <src/util/f77.h>
#include <src/util/taskqueue.h>
#include <src/util/parallel/resources.h>

using namespace std;
using namespace bagel;
using namespace btas;

MP2::MP2(const shared_ptr<const PTree> input, const shared_ptr<const Geometry> g, const shared_ptr<const Reference> ref) : Method(input, g, ref) {

  scf_ = make_shared<RHF>(input, g, ref);
  scf_->compute();
  ref_ = scf_->conv_to_ref();

  cout << endl << "  === DF-MP2 calculation ===" << endl << endl;

  // checks for frozen core
  const bool frozen = idata_->get<bool>("frozen", true);
  ncore_ = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (ncore_) cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;

  if (geom_->df() == nullptr) throw logic_error("MP2 is only implemented with DF");

  // if three is a aux_basis keyword, we use that basis
  abasis_ = to_lower(idata_->get<string>("aux_basis", ""));

}


void MP2::compute() {
  const size_t nbasis = ref_->coeff()->mdim();
  const size_t nocc = ref_->nocc() - ncore_;
  if (nocc < 1)
    throw runtime_error("no correlated electrons");
  if (nbasis <= nocc + ncore_)
    throw runtime_error("no virtuals orbitals");
  const size_t nvirt = nbasis - nocc - ncore_;


  const MatView ocoeff = ref_->coeff()->slice(ncore_, ncore_+nocc);
  const MatView vcoeff = ref_->coeff()->slice(ncore_+nocc, ncore_+nocc+nvirt);

  Timer timer;
  // compute transformed integrals
  shared_ptr<DFDistT> fullt;
  size_t memory_size;

  {
    // first compute half transformed integrals
    shared_ptr<DFHalfDist> half;
    if (abasis_.empty()) {
      half = geom_->df()->compute_half_transform(ocoeff);
      // used later to determine the cache size
      memory_size = half->block(0)->size() * 2;
      mpi__->broadcast(&memory_size, 1, 0);
    } else {
      auto info = make_shared<PTree>(); info->put("df_basis", abasis_);
      auto cgeom = make_shared<Geometry>(*geom_, info, false);
      half = cgeom->df()->compute_half_transform(ocoeff);
      // used later to determine the cache size
      memory_size = cgeom->df()->block(0)->size();
      mpi__->broadcast(&memory_size, 1, 0);
    }

    // second transform for virtual index and rearrange data
    {
      // this is now (naux, nvirt, nocc), distributed by nvirt*nocc. Always naux*nvirt block is localized to one node
      shared_ptr<DFFullDist> full = half->compute_second_transform(vcoeff)->apply_J()->swap();
      auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size(), full->nocc1());
      fullt = make_shared<DFDistT>(full, dist);
    }

    fullt->discard_df();
  }
  assert(fullt->nblocks() == 1);
  const size_t naux = fullt->naux();

  cout << "    * 3-index integral transformation done" << endl;

  // start communication (n fetch behind) - n is determined by memory size
  MP2Cache cache(naux, nocc, nvirt, fullt);

  const int nloop = cache.nloop();
  const int ncache = min(memory_size/(nvirt*nvirt), size_t(20));
  cout << "    * ncache = " << ncache << endl;
  for (int n = 0; n != min(ncache, nloop); ++n)
    cache.block(n, -1);

  // denominator info
  const vector<double> eig(ref_->eig().begin()+ncore_, ref_->eig().end());

  // loop over tasks
  energy_ = 0;
  for (int n = 0; n != nloop; ++n) {
    // take care of data. The communication should be hidden
    if (n+ncache < nloop)
      cache.block(n+ncache, n-1);

    const int i = get<0>(cache.task(n));
    const int j = get<1>(cache.task(n));
    if (i < 0 || j < 0) continue;
    cache.data_wait(n);

    shared_ptr<const Matrix> iblock = cache(i);
    shared_ptr<const Matrix> jblock = cache(j);
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

  // just to double check that all the communition is done
  cache.wait();
  // allreduce energy contributions
  mpi__->allreduce(&energy_, 1);

  cout << "    * assembly done" << endl << endl;
  cout << "      MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << energy_ << setw(10) << setprecision(2) << timer.tick() << endl << endl;

  energy_ += ref_->energy(0);
  cout << "      MP2 total energy:       " << fixed << setw(15) << setprecision(10) << energy_ << endl << endl;
}


