//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ks.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/scf/ks/ks.h>
#include <src/scf/hf/fock.h>
#include <src/prop/multipole.h>
#include <src/util/math/diis.h>

using namespace std;
using namespace bagel;

void KS::compute() {

  Matrix intermediate = *tildex_ % *hcore_ * *tildex_;
  intermediate.diagonalize(eig());
  coeff_ = make_shared<Coeff>(*tildex_ * intermediate);
  shared_ptr<const Matrix> aodensity_ = coeff_->form_density_rhf(nocc_);

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;

  cout << indent << "     - DIIS with orbital gradients will be used." << endl << endl;
  cout << indent << "=== KS iteration (" << name_ << " / " << geom_->basisfile() << ") ===" << endl << indent << endl;

  DIIS<Matrix> diis(diis_size_);

  shared_ptr<Matrix> fock;

  Timer scftime;
  for (int iter = 0; iter != max_iter_; ++iter) {

    // fock operator without DFT xc
    fock = make_shared<Fock<1>>(geom_, hcore_, aodensity_, coeff_->slice(0, nocc_), false /*store*/, true /*rhf*/, func_->scale_ex());

    // add xc
    shared_ptr<const Matrix> xc;
    double exc;
    tie(xc, exc) = grid_->compute_xc(func_, coeff_->slice_copy(0, nocc_));

    energy_ = 0.5*((*hcore_+ *fock) * *aodensity_).trace() + exc + geom_->nuclear_repulsion();

    *fock += *xc;

    auto error_vector = make_shared<const Matrix>(*fock**aodensity_**overlap_ - *overlap_**aodensity_**fock);

    const double error = error_vector->rms();

    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_)
      fock = diis.extrapolate({fock, error_vector});

    {
      auto intermediate = make_shared<Matrix>(*tildex_ % *fock * *tildex_);
      intermediate->diagonalize(eig());
      coeff_ = make_shared<const Coeff>(*tildex_**intermediate);
    }
    aodensity_ = coeff_->form_density_rhf(nocc_);

  }

  // by default we compute dipoles
  if (!geom_->external()) {
    Dipole mu(geom_, aodensity_);
    scf_dipole_ = mu.compute();
  }
}


shared_ptr<const Reference> KS::conv_to_ref() const {
  auto out = make_shared<Reference>(geom_, coeff(), nocc(), 0, coeff()->mdim()-nocc(), vector<double>{energy_});
  out->set_eig(eig_);
  return out;
}

