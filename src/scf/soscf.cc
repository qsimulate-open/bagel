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
}

void SOSCF::initial_guess() {
  shared_ptr<const Matrix> fock = hcore_;
  shared_ptr<Matrix> intermediate = make_shared<Matrix>(*tildex_ % *fock * *tildex_);
  intermediate->diagonalize(eig_.get());
  coeff_ = make_shared<const Coeff> (*tildex_ * *intermediate);
  tie(aodensity_, aodensityD_, aodensityO1_, aodensityO2_) = form_density_soscf();
}

void SOSCF::compute() {

  soeig_ = unique_ptr<double[]> (new double [geom_->nbasis() * 2]);

  initial_guess();

  DIIS<Matrix> diis(diis_size_);

  for (int iter = 0; iter != max_iter_; ++iter) {
    shared_ptr<const Matrix> fock = make_shared<const Fock<1>> (geom_,
                             hcore_, aodensityD_, coeff_->slice(0, nocc_));
    energy_ = 0.5 * ((*hcore_ + *fock) * *aodensityD_).trace() + geom_->nuclear_repulsion();
    auto error_vector = make_shared<const Matrix>(*fock * *aodensityD_ * *overlap_ - 
                                                  *overlap_ * *aodensityD_ * *fock);
    const double error = error_vector->rms();
    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << endl;
    #if 1

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SOSCF iteration converged." << endl << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SOSCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) {
      fock = diis.extrapolate(make_pair(fock, error_vector));
    }

    shared_ptr<Matrix> intermediate = make_shared<Matrix>(*tildex_ % *fock * *tildex_);
    intermediate->diagonalize(eig_.get());
    coeff_ = make_shared<const Coeff>(*tildex_ * *intermediate);
    aodensityD_ = coeff_->form_density_rhf(nocc_);
    #endif

  }
}

tuple<shared_ptr<Matrix>, shared_ptr<Matrix>, shared_ptr<Matrix>, shared_ptr<Matrix>> 
SOSCF::form_density_soscf() const {
  shared_ptr<Matrix> outD = coeff_->form_density_rhf(nocc_);
  // TODO: off-diagonal blocks outO1 and outO2
  shared_ptr<Matrix> outO1 = make_shared<Matrix>(outD->ndim(), outD->mdim());
  shared_ptr<Matrix> outO2 = outO1;
  shared_ptr<Matrix> out = make_shared<Matrix> (2 * outD->ndim(), 2 * outD->mdim());
  out->add_block(1.0, 0, 0, outD->ndim(), outD->mdim(), outD);
  out->add_block(1.0, outD->ndim(), outD->mdim(), outD->ndim(), outD->mdim(), outD);
  out->add_block(1.0, 0, outD->mdim(), outD->ndim(), outD->mdim(), outO1);
  out->add_block(1.0, outD->ndim(), 0, outD->ndim(), outD->mdim(), outO2);
  return make_tuple(out, outD, outO1, outO2);
}

shared_ptr<Matrix> SOSCF::sofock(shared_ptr<const Matrix> fock, 
                   shared_ptr<const Matrix> offd1, shared_ptr<const Matrix> offd2) {
  assert(fock->ndim() == offd1->ndim());
  assert(fock->mdim() == offd1->mdim());
  shared_ptr<Matrix> out = make_shared<Matrix> (2 * fock->ndim(), 2 * fock->mdim());
  out->add_block(1.0, 0, 0, fock->ndim(), fock->mdim(), fock);
  out->add_block(1.0, fock->ndim(), fock->mdim(), fock->ndim(), fock->mdim(), fock);
  out->add_block(1.0, 0, fock->mdim(), fock->ndim(), fock->mdim(), offd1);
  out->add_block(1.0, fock->ndim(), 0, fock->ndim(), fock->mdim(), offd2);
  return out;
}
