//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rohf.cc
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

#include <src/scf/hf/rohf.h>
#include <src/scf/hf/fock.h>
#include <src/prop/multipole.h>
#include <src/util/math/diis.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(ROHF)

void ROHF::compute() {

  initial_guess();

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl;
  cout << endl;
  cout << indent << "    * DIIS with orbital gradients will be used."
            << endl << endl;
  cout << indent << "=== ROHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  // starting SCF iteration

  DIIS<Matrix> diis(diis_size_);
  DIIS<Matrix> diisB(diis_size_);
  Timer scftime;
  for (int iter = 0; iter != max_iter_; ++iter) {

    shared_ptr<const Matrix> fockA = make_shared<const Fock<1>>(geom_, hcore_, aodensity_, coeff_->slice(0,nocc_));
    shared_ptr<const Matrix> fockB = noccB_ ? make_shared<const Fock<1>>(geom_, hcore_, aodensity_, coeffB_->slice(0, noccB_))
                                            : make_shared<const Matrix>(geom_->nbasis(), geom_->nbasis());

    shared_ptr<const Coeff> natorb = get<0>(natural_orbitals());

    auto error_vector = make_shared<Matrix>(*fockA**aodensityA_**overlap_ - *overlap_**aodensityA_**fockA
                                           +*fockB**aodensityB_**overlap_ - *overlap_**aodensityB_**fockB);

    const double error = error_vector->rms();

    energy_ = 0.25*(*aodensityA_*(*hcore_+*fockA)+*aodensityB_*(*hcore_+*fockB)).trace() + geom_->nuclear_repulsion();

    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) {
      shared_ptr<Matrix> tmp_fock = diis.extrapolate({fockA, error_vector});
      shared_ptr<Matrix> intermediateA = make_shared<Matrix>(*natorb % *tmp_fock * *natorb);
                         tmp_fock = diisB.extrapolate({fockB, error_vector});
      shared_ptr<Matrix> intermediateB = make_shared<Matrix>(*natorb % *tmp_fock * *natorb);

      // Specific to ROHF:
      //   here we want to symmetrize closed-virtual blocks
      symmetrize_cv(intermediateA, intermediateB);

      intermediateA->diagonalize(eig());
      intermediateB->diagonalize(eigB());
      coeff_  = make_shared<Coeff>(*natorb * *intermediateA);
      coeffB_ = make_shared<Coeff>(*natorb * *intermediateB);
    } else {
      auto intermediateA = make_shared<Matrix>(*natorb % *fockA * *natorb);
      auto intermediateB = make_shared<Matrix>(*natorb % *fockB * *natorb);

      // Specific to ROHF:
      //   here we want to symmetrize closed-virtual blocks
      symmetrize_cv(intermediateA, intermediateB);

      intermediateA->diagonalize(eig());
      intermediateB->diagonalize(eigB());

      coeff_  = make_shared<Coeff>(*natorb * *intermediateA);
      coeffB_ = make_shared<Coeff>(*natorb * *intermediateB);
    }
    tie(aodensity_, aodensityA_, aodensityB_) = form_density_uhf();

  }

  print_S2("ROHF");
  // by default we compute dipoles
  if (!geom_->external() && multipole_print_) {
    Multipole mu(geom_, aodensity_, multipole_print_);
    scf_dipole_ = mu.compute();
  }
}


void ROHF::symmetrize_cv(shared_ptr<Matrix> fockA, shared_ptr<Matrix> fockB) {
  assert(noccB_ <= nocc_);
  for (int i = 0; i != noccB_; ++i) {
    for (int j = nocc_; j != coeff_->mdim(); ++j) {
      const double dat = (fockA->element(j,i) + fockB->element(j,i)) * 0.5;
      fockA->element(j,i) = fockB->element(j,i) = fockA->element(i,j) = fockB->element(i,j) = dat;
    }
  }
}
