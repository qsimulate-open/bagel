//
// BAGEL - Parallel electron correlation program.
// Filename: mp2.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <stddef.h>
#include <src/mp2/mp2.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <src/util/f77.h>
#include <src/scf/scf.h>
#include <src/smith/prim_op.h>
#include <src/mp2/f12int4.h>
#include <src/util/taskqueue.h>
#include <src/parallel/resources.h>

using namespace std;
using namespace bagel;

MP2::MP2(const multimap<string, string> input, const shared_ptr<const Geometry> g, const shared_ptr<const Reference> ref) : idata_(input), geom_(g) {

  scf_ = std::shared_ptr<SCF<1> >(new SCF<1>(input, g, ref));
  scf_->compute();
  ref_ = scf_->conv_to_ref();

  cout << endl << "  === DF-MP2 calculation ===" << endl << endl;

  // checks for frozen core
  const bool frozen = read_input<bool>(idata_, "frozen", false);
  ncore_ = read_input<int>(idata_, "ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (ncore_) cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;

  ref_->set_ncore(ncore_);

  if (geom_->df() == nullptr) throw logic_error("MP2 is only implemented in DF");

}

namespace bagel {

class MP2AssemTask {
  protected:
    shared_ptr<const DFFullDist> full_;
    const size_t ivirt_;
    const size_t nvirt_;
    const size_t nocc_;
    MP2* mp2_;
  public:
    MP2AssemTask(shared_ptr<const DFFullDist> f, const int iv, const int nv, const int no, MP2* m)
      : full_(f), ivirt_(iv), nvirt_(nv), nocc_(no), mp2_(m) {};

    void compute() {
      shared_ptr<StackMem> stack = resources__->get();
      double* const buf = stack->get(nocc_*nvirt_*nocc_);
      vector<double> eig(mp2_->ref_->eig().begin()+mp2_->ncore_, mp2_->ref_->eig().end());

      // nocc * nvirt * nocc
      unique_ptr<double[]> data = full_->form_4index_1fixed(full_, 1.0, ivirt_);
      copy(data.get(), data.get()+nocc_*nvirt_*nocc_, buf);

      // using SMITH's symmetrizer (src/smith/prim_op.h)
      SMITH::sort_indices<2,1,0,2,1,-1,1>(data.get(), buf, nocc_, nvirt_, nocc_);
      double* tdata = buf;
      for (size_t j = 0; j != nocc_; ++j)
        for (size_t k = 0; k != nvirt_; ++k)
          for (size_t l = 0; l != nocc_; ++l, ++tdata)
            *tdata /= -eig[ivirt_+nocc_]+eig[j]-eig[k+nocc_]+eig[l];

      boost::lock_guard<boost::mutex> lock(mp2_->mut_);
      mp2_->energy_ += ddot_(nocc_*nvirt_*nocc_, data.get(), 1, buf, 1);
      stack->release(nocc_*nvirt_*nocc_, buf);
      resources__->release(stack);
    }
};

}

void MP2::compute() {
  // TODO this factor of 2 is very much error-prone..
  const size_t nocc = geom_->nele()/2 - ncore_;
  if (nocc < 1) throw runtime_error("no correlated electrons");
  const size_t nvirt = geom_->nbasis() - nocc - ncore_;
  if (nvirt < 1) throw runtime_error("no virtuals orbitals");
  assert(geom_->nbasis() == ref_->coeff()->mdim());

  const size_t nbasis = geom_->nbasis();

  const double* const coeff = ref_->coeff()->data() + ncore_*nbasis;
  const double* const vcoeff = coeff + nocc*nbasis;

  auto tp1 = chrono::high_resolution_clock::now();

  // first compute half transformed integrals
  shared_ptr<DFHalfDist> half = geom_->df()->compute_half_transform(coeff, nocc);
  // second transform for virtual index
  // this is now (naux, nocc, nvirt)
  shared_ptr<DFFullDist> full = half->compute_second_transform(vcoeff, nvirt)->apply_J();

  cout << "    * 3-index integral transformation done" << endl;

  // assemble
  energy_ = 0.0;
  vector<MP2AssemTask> task;
  for (size_t i = 0; i != nvirt; ++i) {
    MP2AssemTask tmp(full, i, nvirt, nocc, this);
    task.push_back(tmp);
  }

  TaskQueue<MP2AssemTask> tq(task);
  tq.compute(resources__->max_num_threads());

  auto tp2 = chrono::high_resolution_clock::now();
  cout << "    * assembly done" << endl << endl;
  cout << "      MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << energy_
                      << setw(10) << setprecision(2) << chrono::duration_cast<chrono::milliseconds>(tp2-tp1).count()*0.001 << endl << endl;

  energy_ += ref_->energy();
  cout << "      MP2 total energy:       " << fixed << setw(15) << setprecision(10) << energy_ << endl << endl;

  // check if F12 is requested.
  const bool do_f12 = read_input<bool>(idata_, "f12", false);
  if (do_f12) {
    const double gamma = read_input<double>(idata_, "gamma", 1.5);
    cout << "    * F12 calculation requested with gamma = " << setprecision(2) << gamma << endl;
#if 0
    shared_ptr<F12Int> f12int(new F12Int(idata_, geom_, ref_, gamma, ncore_));
#else
    shared_ptr<F12Ref> f12ref(new F12Ref(geom_, ref_, ncore_, gamma));
    f12ref->compute();
#endif

  }
}


