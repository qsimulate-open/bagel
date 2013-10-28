//
// BAGEL - Parallel electron correlation program.
// Filename: soscf.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu> 
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

#include <src/scf/soscf.h>
#include <src/math/diis.h>
#include <src/prop/multipole.h>
#include <src/scf/atomicdensities.h>


using namespace std;
using namespace bagel;

SOSCF::SOSCF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
 : SCF_base(idata, geom, re) {
  cout << indent << "*** Two-component ECP-SCF ***" << endl << endl;
  soeig_ = unique_ptr<double[]> (new double[2 * geom_->nbasis()]);
  sohcore_ = make_shared<SOhcore>(geom_, hcore_);
}

void SOSCF::initial_guess() {
  sooverlap_ = SOSCF::sooverlap();
  sotildex_ = SOSCF::sotildex();

  shared_ptr<const Matrix> fock = hcore_;
#if 1
  shared_ptr<Matrix> intermediate = make_shared<Matrix>(*tildex_ % *fock * *tildex_);
  intermediate->diagonalize(eig_.get());
  coeff_ = make_shared<const Coeff> (*tildex_ * *intermediate);
  aodensityD_ = coeff_->form_density_rhf(nocc_);
#endif

#if 1
  shared_ptr<const Matrix> sofock = SOSCF::sofock(fock, sohcore_);
  shared_ptr<Matrix> sointermediate = make_shared<Matrix>(*sotildex_ % *sofock * *sotildex_);
  sointermediate->diagonalize(soeig_.get());
  socoeff_ = make_shared<Coeff>(*sotildex_ * *sointermediate);
  aodensity_ = socoeff_->form_density_rhf(2 * nocc_);
#endif
}

void SOSCF::compute() {
#if 1
  initial_guess();

  DIIS<Matrix> diis(diis_size_);

for (int iter = 0; iter != max_iter_; ++iter) {
    shared_ptr<const Matrix> fock = make_shared<const Fock<1>> (geom_,
                             hcore_, aodensityD_, coeff_->slice(0, nocc_));
    #if 1
    shared_ptr<const Matrix> sofock = SOSCF::sofock(fock, sohcore_);
    energy_ = 0.25 * ((*sohcore_ + *sofock) * *aodensity_).trace() + geom_->nuclear_repulsion();
    auto error_vector = make_shared<const Matrix>(*sofock * *aodensity_ * *sooverlap_ - 
                                                  *sooverlap_ * *aodensity_ * *sofock);
    const double error = error_vector->rms();
    #endif
    #if 0
    energy_ = 0.5 * ((*hcore_ + *fock) * *aodensityD_).trace() + geom_->nuclear_repulsion();
    auto error_vector = make_shared<const Matrix>(*fock * *aodensityD_ * *overlap_ - 
                                                  *overlap_ * *aodensityD_ * *fock);
    const double error = error_vector->rms();
    #endif

    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << endl;
    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SOSCF iteration converged." << endl << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SOSCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) {
      sofock = diis.extrapolate(make_pair(sofock, error_vector));
    }

    #if 1
    shared_ptr<Matrix> sointermediate = make_shared<Matrix>(*sotildex_ % *sofock * *sotildex_);
    sointermediate->diagonalize(soeig_.get());
    shared_ptr<Coeff> socoeff_ = make_shared<Coeff>(*sotildex_ * *sointermediate);
    aodensity_ = socoeff_->form_density_rhf(2 * nocc_);
    #endif

    #if 1
    shared_ptr<Matrix> intermediate = make_shared<Matrix>(*tildex_ % *fock * *tildex_);
    intermediate->diagonalize(eig_.get());
    coeff_ = make_shared<Coeff>(*tildex_ * *intermediate);
    aodensityD_ = coeff_->form_density_rhf(nocc_);
    #endif
  }
#endif
}

shared_ptr<Matrix> SOSCF::sofock(shared_ptr<const Matrix> fock, shared_ptr<const Matrix> sohcore) {
  assert(sohcore->ndim() == 2 * fock->ndim());
  assert(sohcore->mdim() == 2 * fock->mdim());
  shared_ptr<Matrix> out = make_shared<Matrix>(sohcore->ndim(), sohcore->mdim());
  out->add_block(1.0, 0, 0, fock->ndim(), fock->mdim(), fock);
  out->add_block(1.0, fock->ndim(), fock->mdim(), fock->ndim(), fock->mdim(), fock);
  return out;
}

shared_ptr<Matrix> SOSCF::sotildex() {
  shared_ptr<Matrix> out = make_shared<Matrix>(2 * tildex_->ndim(), 2 * tildex_->mdim());
  out->zero();
  out->add_block(1.0, 0, 0, tildex_->ndim(), tildex_->mdim(), tildex_);
  out->add_block(1.0, tildex_->ndim(), tildex_->mdim(), tildex_->ndim(), tildex_->mdim(), tildex_);
  return out;
}

shared_ptr<Matrix> SOSCF::sooverlap() {
  shared_ptr<Matrix> out = make_shared<Matrix>(2 * overlap_->ndim(), 2 * overlap_->mdim());
  out->zero();
  out->add_block(1.0, 0, 0, overlap_->ndim(), overlap_->mdim(), overlap_);
  out->add_block(1.0, overlap_->ndim(), overlap_->mdim(), overlap_->ndim(), overlap_->mdim(), overlap_);
  return out;
}

