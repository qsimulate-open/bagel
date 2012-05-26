//
// Newint - Parallel electron correlation program.
// Filename: mp2.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


// implements the MP2-F12 theory

#include <src/mp2/mp2.h>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/scf/scf.h>
#include <src/smith/prim_op.h>
#include <src/mp2/f12int4.h>

using namespace std;

MP2::MP2(const multimap<string, string> input, const shared_ptr<const Geometry> g) : idata_(input), geom_(g) {

  scf_ = std::shared_ptr<SCF<1> >(new SCF<1>(input, g));
  scf_->compute();
  ref_ = scf_->conv_to_ref();

  cout << endl << "  === DF-MP2 calculation ===" << endl << endl;

  // checks for frozen core
  const bool frozen = read_input<bool>(idata_, "frozen", false);
  ncore_ = read_input<int>(idata_, "ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (ncore_) cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;

  ref_->set_ncore(ncore_);

  if (!geom_->df()) throw logic_error("MP2 is only implemented in DF");

}

void MP2::compute() {
  // TODO this factor of 2 is very much error-prone..
  const size_t nocc = geom_->nele() / 2 - ncore_;
  if (nocc < 1) throw runtime_error("no correlated electrons"); 
  const size_t nvirt = geom_->nbasis() - nocc - ncore_;
  if (nvirt < 1) throw runtime_error("no virtuals orbitals"); 
  assert(geom_->nbasis() == ref_->coeff()->mdim());

  const size_t nbasis = geom_->nbasis();

  const double* const coeff = ref_->coeff()->data() + ncore_*nbasis;
  const double* const vcoeff = coeff + nocc*nbasis;

  const long time = ::clock();

  // first compute half transformed integrals
  shared_ptr<DF_Half> half = geom_->df()->compute_half_transform(coeff, nocc);  
  // second transform for virtual index
  // this is now (naux, nocc, nvirt)
  shared_ptr<DF_Full> full = half->compute_second_transform(vcoeff, nvirt)->apply_J();

  cout << "    * 3-index integral transformation done" << endl;

  // assemble
  unique_ptr<double[]> buf(new double[nocc*nvirt*nocc]);
  vector<double> eig_tm = ref_->eig();
  vector<double> eig(eig_tm.begin()+ncore_, eig_tm.end());

  // TODO in priciple this should run over occupied (for optimal implementations)...

  double sum = 0.0;
  for (size_t i = 0; i != nvirt; ++i) {
    // nocc * nvirt * nocc
    unique_ptr<double[]> data = full->form_4index(full, i); 
    copy(data.get(), data.get()+nocc*nvirt*nocc, buf.get());

    // using SMITH's symmetrizer (src/smith/prim_op.h)
    SMITH::sort_indices<2,1,0,2,1,-1,1>(data, buf, nocc, nvirt, nocc);
    double* tdata = buf.get();
    for (size_t j = 0; j != nocc; ++j) {
      for (size_t k = 0; k != nvirt; ++k) {
        for (size_t l = 0; l != nocc; ++l, ++tdata) {
          *tdata /= -eig[i+nocc]+eig[j]-eig[k+nocc]+eig[l]; 
        }
      }
    }
    sum += ddot_(nocc*nvirt*nocc, data, 1, buf, 1);
  }

  const double elapsed = (::clock()-time)/static_cast<double>(CLOCKS_PER_SEC); 
  cout << "    * assembly done" << endl << endl;
  cout << "      MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << sum
                                                    << setw(10) << setprecision(2) << elapsed << endl << endl;

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


