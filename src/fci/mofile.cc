//
// Newint - Parallel electron correlation program.
// Filename: mofile.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <src/util/f77.h>
#include <src/fci/mofile.h>
#include <src/scf/scf.h>

using namespace std;

typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<Shell> RefShell;

MOFile::MOFile(const shared_ptr<const Geometry> geom, const shared_ptr<Reference> ref) : geom_(geom), ref_(ref),
    core_fock_(new double[geom->nbasis()*geom->nbasis()]) {

  do_df_ = geom->df().get();
  if (!do_df_) throw runtime_error("for the time being I gave up maintaining non-DF codes.");

}

MOFile::~MOFile() {
}



double MOFile::create_Jiiii(const int nstart, const int nfence, const int ncore) {

  // first compute all the AO integrals in core
  nocc_ = nfence - nstart;
  nbasis_ = geom_->nbasis(); // size_t quantity :-)
  const int nbasis = nbasis_;
  const int nocc = nocc_;

  // some stuffs for blas
  double* cdata = ref_->coeff()->data() + nstart*nbasis;
  const int mm = nocc*nocc;

  // one electron part
  double core_energy = 0.0;
  {
    unique_ptr<double[]> aobuff(new double[nbasis*nbasis]);
    shared_ptr<Fock<1> > fock0(new Fock<1>(geom_, ref_->hcore()));
    if (nstart != 0) {
      shared_ptr<Matrix1e> den = ref_->coeff()->form_density_rhf(ncore);
      fock0 = shared_ptr<Fock<1> >(new Fock<1>(geom_, fock0, den, ref_->schwarz()));
      core_energy = (*den * (*ref_->hcore()+*fock0)).trace() * 0.5;
      dcopy_(nbasis*nbasis, fock0->data(), 1, core_fock_ptr(), 1);
    }
    fock0->fill_upper();
    dgemm_("n","n",nbasis,nocc,nbasis,1.0,fock0->data(),nbasis,cdata,nbasis,0.0,aobuff.get(),nbasis);

    mo1e_ = unique_ptr<double[]>(new double[nocc*nocc]);
    dgemm_("t","n",nocc,nocc,nbasis,1.0,cdata,nbasis,aobuff.get(),nbasis,0.0,mo1e_ptr(),nocc);
  }


  //
  // two electron part.
  //

  shared_ptr<const DensityFit> dff = geom_->df();

  // first half transformation
  shared_ptr<DF_Half> half = dff->compute_half_transform(cdata, nocc);

  // second index transformation and (D|ii) = J^-1/2_DE (E|ii)
  shared_ptr<DF_Full> buf = half->compute_second_transform(cdata, nocc)->apply_J();

  // assembles (ii|ii) = (ii|D)(D|ii)
  unique_ptr<double[]> buf2e(new double[nocc*nocc*nocc*nocc]);
  buf->form_4index(buf2e);

  // we want to store half-transformed quantity for latter convenience
  mo2e_1ext_size_ = nocc*dff->naux()*nbasis;
  mo2e_1ext_ = half;


  // mo2e is compressed
  sizeij_ = nocc*(nocc+1)/2;
  mo2e_ = unique_ptr<double[]>(new double[sizeij_*sizeij_]);

  int ijkl = 0;
  for (int i = 0; i != nocc; ++i) {
    for (int j = 0; j <= i; ++j) {
      const int ijo = (j + i*nocc)*nocc*nocc;
      for (int k = 0; k != nocc; ++k) {
        for (int l = 0; l <= k; ++l, ++ijkl) {
          mo2e_[ijkl] = buf2e[l+k*nocc+ijo]; 
        }
      }
    }
  }

  // h'kl = hkl - 0.5 sum_j (kj|jl)
  unique_ptr<double[]> buf3(new double[sizeij_]);
  int ij = 0;
  for (int i=0; i!=nocc; ++i) {
    for (int j=0; j<=i; ++j, ++ij) {
      buf3[ij] = mo1e_[j+i*nocc];
      for (int k=0; k!=nocc; ++k) {
        buf3[ij] -= 0.5*buf2e[(k+j*nocc)*mm+(k+i*nocc)];
      }
    }
  }
  mo1e_ = move(buf3);
  return core_energy;
}


void MOFile::update_1ext_ints(const vector<double>& coeff) {
  // in the case of no DF
  shared_ptr<DF_Half> buf = mo2e_1ext_->clone();

  // half transformed DF is rotated.
  const int naux = geom_->df()->naux();
  for (int i = 0; i != nbasis_; ++i) 
    dgemm_("N", "N", naux, nocc_, nocc_, 1.0, mo2e_1ext_->data()+i*naux*nocc_, naux, &(coeff[0]), nocc_, 0.0, buf->data()+i*naux*nocc_, naux); 
  mo2e_1ext_ = buf;

}
