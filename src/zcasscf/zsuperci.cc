//
// BAGEL - Parallel electron correlation program.
// Filename: zuperci.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Jefferson Bates <jefferson.bates@northwestern.edu>
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


#include <src/zcasscf/zqvec.h>
#include <src/math/step_restrict_bfgs.h>
#include <src/rel/dfock.h>
#include <src/zcasscf/zsuperci.h>
#include <src/rel/reloverlap.h>

using namespace std;
using namespace bagel;

void ZSuperCI::compute() {
    // TODO : DIIS needs to be introduced eventually

  // ============================
  // macro iteration from here
  // ============================
   Timer timer;

  // intialize coefficients
  init_kramers_coeff();

  cout << setprecision(8) << " kramers restricted re part rms = " << coeff_->get_real_part()->rms() << endl;
  cout << setprecision(8) << " kramers restricted im part rms = " << coeff_->get_imag_part()->rms() << endl;

  if (nact_)
    fci_->update(coeff_);

  for (int iter = 0; iter != 1; ++iter) {

    // first perform CASCI to obtain RDMs
    if (nact_) {
      mute_stdcout(/*fci*/true);
      if (iter) fci_->update(coeff_);
      cout << " Executing FCI calculation in Cycle " << iter << endl;
      fci_->compute();
      cout << " Computing RDMs from FCI calculation " << endl;
      fci_->compute_rdm12();
      resume_stdcout();
    }

  }
}


pair<shared_ptr<ZMatrix>, vector<double>> ZSuperCI::make_natural_orbitals(shared_ptr<const ZMatrix> rdm1) const {
  // input should be 1rdm in kramers format
  shared_ptr<ZMatrix> tmp = rdm1->copy();

  unique_ptr<double[]> vec(new double[rdm1->ndim()]);
  zquatev_(tmp->ndim(), tmp->data(), vec.get());

  map<int,int> emap;
  auto buf2 = tmp->clone();
  vector<double> vec2(tmp->ndim());
  // sort eigenvectors so that buf is close to a unit matrix
  // target column
  for (int i = 0; i != tmp->ndim(); ++i) {
    // first find the source column
    tuple<int, double> max = make_tuple(-1, 0.0);
    for (int j = 0; j != tmp->ndim(); ++j)
      if (sqrt(real(tmp->element(i,j)*tmp->get_conjg()->element(i,j))) > get<1>(max))
        max = make_tuple(j, sqrt(real(tmp->element(i,j)*tmp->get_conjg()->element(i,j))));

    // register to emap
    if (emap.find(get<0>(max)) != emap.end()) throw logic_error("this should not happen. make_natural_orbitals()");
    emap.insert(make_pair(get<0>(max), i));

    // copy to the target
    copy_n(tmp->element_ptr(0,get<0>(max)), tmp->ndim(), buf2->element_ptr(0,i));
    vec2[i] = vec[get<0>(max)];
  }

  // fix the phase
  for (int i = 0; i != tmp->ndim(); ++i) {
    if (real(buf2->element(i,i)) < 0.0)
      blas::scale_n(-1.0, buf2->element_ptr(0,i), tmp->ndim());
  }
  // copy eigenvalues TODO: change to blas
  for (int i=0; i!=tmp->ndim()/2; ++i)
    vec2[tmp->ndim()/2 + i] = vec2[i];

  return make_pair(buf2, vec2);
}


shared_ptr<const ZMatrix> ZSuperCI::natorb_rdm1_transform(const shared_ptr<ZMatrix> coeff, shared_ptr<const ZMatrix> rdm1) {
  shared_ptr<ZMatrix> tmp = rdm1->clone();
  const complex<double>* start = coeff->data();
  int ndim = coeff->ndim();
  unique_ptr<complex<double>[]> buf(new complex<double>[ndim*ndim]);
  zgemm3m_("N", "N", ndim, ndim, ndim, 1.0, rdm1->data(), ndim, start, ndim, 0.0, buf.get(), ndim);
  zgemm3m_("T", "N", ndim, ndim, ndim, 1.0, start, ndim, buf.get(), ndim, 0.0, tmp->data(), ndim);
  auto out = make_shared<const ZMatrix>(*tmp); 
  return out;
}
