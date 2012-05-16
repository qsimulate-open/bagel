//
// Newint - Parallel electron correlation program.
// Filename: rohf.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/scf/rohf.h>

using namespace std;

void ROHF::compute() {

  string indent = "  ";
  shared_ptr<Fock<1> > hcore_fock;
  {
    hcore_fock = shared_ptr<Fock<1> >(new Fock<1>(geom_, hcore_));
   
    Matrix1e intermediate = *tildex_ % *hcore_fock * *tildex_;
    intermediate.diagonalize(eig());
    coeff_ = shared_ptr<Coeff>(new Coeff(*tildex_ * intermediate));
    coeffB_ = shared_ptr<Coeff>(new Coeff(*coeff_)); // since this is obtained with hcore
    tie(aodensity_, aodensityA_, aodensityB_) = form_density_uhf();
  }

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl;
  cout << endl; 
  cout << indent << "    * DIIS with " << (density_change_ ? "density changes" : "orbital gradients") << " will be used."
            << endl << endl;
  cout << indent << "=== ROHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  // starting SCF iteration
  eigB_ = unique_ptr<double[]>(new double[coeff_->mdim()]);

  DIIS<Matrix1e> diis(5);
  DIIS<Matrix1e> diisB(5);
  for (int iter = 0; iter != max_iter_; ++iter) {
    int start = ::clock();

    shared_ptr<Fock<1> > fockA(new Fock<1>(geom_, hcore_fock, aodensity_, aodensityA_, schwarz_));
    shared_ptr<Fock<1> > fockB(new Fock<1>(geom_, hcore_fock, aodensity_, aodensityB_, schwarz_));

    shared_ptr<const Coeff> natorb = get<0>(natural_orbitals());

    shared_ptr<Matrix1e> intermediateA(new Matrix1e(*natorb % *fockA * *natorb));
    shared_ptr<Matrix1e> intermediateB(new Matrix1e(*natorb % *fockB * *natorb));

    // Specific to ROHF:
    //   here we want to symmetrize closed-virtual blocks
    symmetrize_cv(intermediateA, intermediateB);

    intermediateA->diagonalize(eig());
    intermediateB->diagonalize(eigB());

    shared_ptr<Matrix1e> new_density, new_densityA, new_densityB;
    if (density_change_)
      tie(new_density, new_densityA, new_densityB) = form_density_uhf();

    shared_ptr<Matrix1e> error_vector(new Matrix1e(
      density_change_ ? (*new_density - *aodensity_) : (*fockA**aodensityA_**overlap_ - *overlap_**aodensityA_**fockA
                                                       +*fockB**aodensityB_**overlap_ - *overlap_**aodensityB_**fockB)));
    
    const double error = error_vector->rms();

    energy_ = 0.5*(*aodensity_ * *hcore_).trace() + geom_->nuclear_repulsion();
    for (int i = 0; i != this->nocc(); ++i)  energy_ += eig_[i]  * 0.5;
    for (int i = 0; i != this->noccB(); ++i) energy_ += eigB_[i] * 0.5;

    int end = ::clock();
    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2)
                                      << (end - start) / static_cast<double>(CLOCKS_PER_SEC) << endl; 

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) {
      shared_ptr<Matrix1e> tmp_fock = diis.extrapolate(make_pair(fockA, error_vector));
      shared_ptr<Matrix1e> intermediateA(new Matrix1e(*natorb % *tmp_fock * *natorb));
                           tmp_fock = diisB.extrapolate(make_pair(fockB, error_vector));
      shared_ptr<Matrix1e> intermediateB(new Matrix1e(*natorb % *tmp_fock * *natorb));

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
  }

  print_S2("ROHF");
}


void ROHF::symmetrize_cv(shared_ptr<Matrix1e> fockA, shared_ptr<Matrix1e> fockB) {
  assert(noccB_ <= nocc_);
  for (int i = 0; i != noccB_; ++i) {
    for (int j = nocc_; j != geom_->nbasis(); ++j) {
      const double dat = (fockA->element(j,i) + fockB->element(j,i)) * 0.5;
      fockA->element(j,i) = fockB->element(j,i) = fockA->element(i,j) = fockB->element(i,j) = dat;
    }
  }
}
