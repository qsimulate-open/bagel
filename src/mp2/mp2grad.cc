//
// Newint - Parallel electron correlation program.
// Filename: mp2grad.cc
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


#include <src/mp2/mp2grad.h>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/smith/prim_op.h>
#include <src/prop/dipole.h>

using namespace std;

MP2Grad::MP2Grad(const multimap<string, string> input, const shared_ptr<Geometry> g, shared_ptr<Reference> r) : MP2(input, g, r) {

}

void MP2Grad::compute() {
  const long time = ::clock();

  // since this is only for closed shell
  const size_t naux = geom_->naux();
  const size_t nocc = geom_->nele() / 2 - ncore_;
  if (nocc < 1) throw runtime_error("no correlated electrons"); 
  const size_t nvirt = geom_->nbasis() - nocc - ncore_;
  if (nvirt < 1) throw runtime_error("no virtuals orbitals"); 
  assert(geom_->nbasis() == ref_->coeff()->mdim());

  const size_t nbasis = geom_->nbasis();

  const double* const coeff = ref_->coeff()->data() + ncore_*nbasis;
  const double* const vcoeff = coeff + nocc*nbasis;

  // first compute half transformed integrals
  shared_ptr<DF_Half> half = geom_->df()->compute_half_transform(coeff, nocc);  
  // second transform for virtual index
  // this is now (naux, nocc, nvirt)
  shared_ptr<DF_Full> full = half->compute_second_transform(vcoeff, nvirt)->apply_J();
  shared_ptr<DF_Full> bv = full->apply_J();
  shared_ptr<DF_Full> gia = bv->clone();

  cout << "    * 3-index integral transformation done" << endl;

  // assemble
  unique_ptr<double[]> buf(new double[nocc*nvirt*nocc]); // it is implicitly assumed that o^2v can be kept in core in each node
  unique_ptr<double[]> buf2(new double[nocc*nvirt*nocc]);
  vector<double> eig_tm = ref_->eig();
  vector<double> eig(eig_tm.begin()+ncore_, eig_tm.end());

  shared_ptr<Matrix1e> dmp2(new Matrix1e(geom_));
  double* optr = dmp2->element_ptr(ncore_, ncore_);
  double* vptr = dmp2->element_ptr(ncore_+nocc, ncore_+nocc);

  double sum = 0.0;
  for (size_t i = 0; i != nvirt; ++i) {
    // nocc * nvirt * nocc
    unique_ptr<double[]> data = full->form_4index(full, i); 
    copy(data.get(), data.get()+nocc*nvirt*nocc, buf.get());
    copy(data.get(), data.get()+nocc*nvirt*nocc, buf2.get());

    // using SMITH's symmetrizer (src/smith/prim_op.h)
    SMITH::sort_indices<2,1,0,2,1,-1,1>(data, buf, nocc, nvirt, nocc);
    double* tdata = buf.get();
    double* bdata = buf2.get();
    for (size_t j = 0; j != nocc; ++j) {
      for (size_t k = 0; k != nvirt; ++k) {
        for (size_t l = 0; l != nocc; ++l, ++tdata, ++bdata) {
          const double denom = 1.0 / (-eig[i+nocc]+eig[j]-eig[k+nocc]+eig[l]);
          *tdata *= denom;
          *bdata *= denom; 
        }
      }
    }
    sum += ddot_(nocc*nvirt*nocc, data, 1, buf, 1);

    // form Gia : TODO distribute
    // Gia(D|ic) = BV(D|ja) G_c(ja|i)
    dgemm_("N", "N", naux, nocc, nvirt*nocc, 1.0, bv->data(), naux, buf.get(), nocc*nvirt, 0.0, gia->data()+i*nocc*naux, naux); 

    // G(ja|ic) -> G_c(a,ij) 
    SMITH::sort_indices<1,2,0,0,1,1,1>(buf, data, nocc, nvirt, nocc);
    // T(jb|ic) -> T_c(b,ij) 
    SMITH::sort_indices<1,2,0,0,1,1,1>(buf2, buf, nocc, nvirt, nocc);
    // D_ab = G(ja|ic) T(jb|ic)
    dgemm_("N", "T", nvirt, nvirt, nocc*nocc, 2.0, buf.get(), nvirt, data.get(), nvirt, 1.0, vptr, nbasis); 
    // D_ij = - G(ja|kc) T(ia|kc)
    dgemm_("T", "N", nocc, nocc, nvirt*nocc, -2.0, buf.get(), nvirt*nocc, data.get(), nvirt*nocc, 1.0, optr, nbasis); 
  }

  // L''aq = 2 Gia(D|ia) (D|iq)
  unique_ptr<double[]> laq = gia->form_2index(half, 2.0);
  unique_ptr<double[]> lai(new double[nocc*nvirt]);
  dgemm_("N", "N", nvirt, nocc, nbasis, 1.0, laq.get(), nvirt, coeff, nbasis, 0.0, lai.get(), nvirt); 

  // Gip = Gia(D|ia) C+_ap
  shared_ptr<DF_Half> gip = gia->back_transform(vcoeff);
  // Liq = 2 Gip(D) * (D|pq)
  unique_ptr<double[]> lip = gip->form_2index(geom_->df(), 2.0);
  unique_ptr<double[]> lia(new double[nocc*nvirt]);
  dgemm_("N", "N", nocc, nvirt, nbasis, 1.0, lip.get(), nocc, vcoeff, nbasis, 0.0, lia.get(), nocc); 

  // printout right hand side
  cout << "  -- printing out the right hand side --" << endl;
  for (int a = 0; a != nvirt; ++a) {
    for (int i = 0; i != nocc; ++i) {
      cout << setw(15) << setprecision(10) << lai[a+nvirt*i] - lia[i+nocc*a];
    }
    cout << endl;
  }
  cout << "  --------------------------------------" << endl;

  const double elapsed = (::clock()-time)/static_cast<double>(CLOCKS_PER_SEC); 
  cout << "    * assembly done" << endl << endl;
  cout << "      MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << sum
                                                    << setw(10) << setprecision(2) << elapsed << endl << endl;

  // computes dipole mements
  shared_ptr<Matrix1e> dmp2ao(new Matrix1e(*ref_->coeff() * *dmp2 ^ *ref_->coeff()));
  Dipole dipole(geom_, dmp2ao);
  dipole.compute();

}


