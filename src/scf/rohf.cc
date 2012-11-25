//
// BAGEL - Parallel electron correlation program.
// Filename: rohf.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <chrono>
#include <src/scf/rohf.h>
#include <src/prop/dipole.h>
#include <src/parallel/paramatrix.h>

using namespace std;
using namespace std::chrono;
using namespace bagel;

void ROHF::compute() {

  string indent = "  ";
  shared_ptr<Fock<1> > hcore_fock(new Fock<1>(geom_, hcore_));

  if (coeff_ == nullptr || coeffB_ == nullptr) {
    ParaMatrix intermediate = *tildex_ % *hcore_fock * *tildex_;
    intermediate.diagonalize(eig());
    coeff_ = shared_ptr<Coeff>(new Coeff(*tildex_ * intermediate));
    coeffB_ = shared_ptr<Coeff>(new Coeff(*coeff_)); // since this is obtained with hcore
  }
  tie(aodensity_, aodensityA_, aodensityB_) = form_density_uhf();

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl;
  cout << endl;
  cout << indent << "    * DIIS with " << (density_change_ ? "density changes" : "orbital gradients") << " will be used."
            << endl << endl;
  cout << indent << "=== ROHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  // starting SCF iteration
  eigB_ = unique_ptr<double[]>(new double[coeff_->mdim()]);

  DIIS<Matrix> diis(5);
  DIIS<Matrix> diisB(5);
  for (int iter = 0; iter != max_iter_; ++iter) {
    auto tp1 = high_resolution_clock::now();

    shared_ptr<Fock<1> > fockA(new Fock<1>(geom_, hcore_fock, aodensity_, aodensityA_, schwarz_));
    shared_ptr<Fock<1> > fockB(new Fock<1>(geom_, hcore_fock, aodensity_, aodensityB_, schwarz_));

    shared_ptr<const Coeff> natorb = get<0>(natural_orbitals());

    shared_ptr<ParaMatrix> intermediateA(new ParaMatrix(*natorb % *fockA * *natorb));
    shared_ptr<ParaMatrix> intermediateB(new ParaMatrix(*natorb % *fockB * *natorb));

    // Specific to ROHF:
    //   here we want to symmetrize closed-virtual blocks
    symmetrize_cv(intermediateA, intermediateB);

    intermediateA->diagonalize(eig());
    intermediateB->diagonalize(eigB());

    shared_ptr<Matrix> new_density, new_densityA, new_densityB;
    if (density_change_)
      tie(new_density, new_densityA, new_densityB) = form_density_uhf();

    shared_ptr<Matrix> error_vector(new Matrix(
      density_change_ ? (*new_density - *aodensity_) : (*fockA**aodensityA_**overlap_ - *overlap_**aodensityA_**fockA
                                                       +*fockB**aodensityB_**overlap_ - *overlap_**aodensityB_**fockB)));

    const double error = error_vector->rms();

    energy_ = 0.5*(*aodensity_ * *hcore_).trace() + geom_->nuclear_repulsion();
    for (int i = 0; i != this->nocc(); ++i)  energy_ += eig_[i]  * 0.5;
    for (int i = 0; i != this->noccB(); ++i) energy_ += eigB_[i] * 0.5;

    auto tp2 = high_resolution_clock::now();
    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2)
                                      << duration_cast<milliseconds>(tp2-tp1).count()*0.001 << endl;

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) {
      shared_ptr<Matrix> tmp_fock = diis.extrapolate(make_pair(fockA, error_vector));
      shared_ptr<ParaMatrix> intermediateA(new ParaMatrix(*natorb % *tmp_fock * *natorb));
                         tmp_fock = diisB.extrapolate(make_pair(fockB, error_vector));
      shared_ptr<ParaMatrix> intermediateB(new ParaMatrix(*natorb % *tmp_fock * *natorb));

      // Specific to ROHF:
      //   here we want to symmetrize closed-virtual blocks
      symmetrize_cv(intermediateA, intermediateB);

      intermediateA->diagonalize(eig());
      intermediateB->diagonalize(eigB());
      coeff_  = shared_ptr<Coeff>(new Coeff(*natorb**intermediateA));
      coeffB_ = shared_ptr<Coeff>(new Coeff(*natorb**intermediateB));
    } else {
      coeff_  = shared_ptr<Coeff>(new Coeff(*natorb * *intermediateA));
      coeffB_ = shared_ptr<Coeff>(new Coeff(*natorb * *intermediateB));
    }
    tie(aodensity_, aodensityA_, aodensityB_) = form_density_uhf();

    // need to make all the node consistent
    mpi__->broadcast(aodensityA_->data(), aodensityA_->size(), 0);
    mpi__->broadcast(aodensityB_->data(), aodensityB_->size(), 0);
    mpi__->broadcast(aodensity_->data(), aodensity_->size(), 0);
  }

  print_S2("ROHF");
  // by default we compute dipoles
  if (!geom_->external()) {
    Dipole mu(geom_, aodensity_);
    mu.compute();
  }
}


void ROHF::symmetrize_cv(shared_ptr<Matrix> fockA, shared_ptr<Matrix> fockB) {
  assert(noccB_ <= nocc_);
  for (int i = 0; i != noccB_; ++i) {
    for (int j = nocc_; j != geom_->nbasis(); ++j) {
      const double dat = (fockA->element(j,i) + fockB->element(j,i)) * 0.5;
      fockA->element(j,i) = fockB->element(j,i) = fockA->element(i,j) = fockB->element(i,j) = dat;
    }
  }
}
