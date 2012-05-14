//
// Newint - Parallel electron correlation program.
// Filename: uhf.cc
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


#include <src/scf/uhf.h>

using namespace std;

void UHF::compute() {

  string indent = "  ";
  shared_ptr<Fock<1> > hcore_fock;
  {
    hcore_fock = shared_ptr<Fock<1> >(new Fock<1>(geom_, hcore_));
   
    Matrix1e intermediate = *tildex_ % *hcore_fock * *tildex_;
    intermediate.diagonalize(eig());
    coeff_ = shared_ptr<Coeff>(new Coeff(*tildex_ * intermediate));
    coeffB_ = shared_ptr<Coeff>(new Coeff(*coeff_));
    aodensity_ = form_density_rhf();
  }

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl;
  cout << endl; 
  cout << indent << "    * DIIS with " << (density_change_ ? "density changes" : "orbital gradients") << " will be used."
            << endl << endl;
  cout << indent << "=== UHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  // starting SCF iteration

  DIIS<Matrix1e> diis(5);
  DIIS<Matrix1e> diisB(5);
  for (int iter = 0; iter != max_iter_; ++iter) {
    int start = ::clock();

    shared_ptr<Fock<1> > fock(new Fock<1>(geom_, hcore_fock, aodensity_, aodensity_, schwarz_));

    Matrix1e intermediate = *coeff_ % *fock * *coeff_;

    intermediate.diagonalize(eig());
    coeff_ = shared_ptr<Coeff>(new Coeff((*coeff_) * intermediate));
    shared_ptr<Matrix1e> new_density = form_density_rhf();

    shared_ptr<Matrix1e> error_vector(new Matrix1e(
      density_change_ ? (*new_density - *aodensity_) : (*fock**aodensity_**overlap_ - *overlap_**aodensity_**fock)
    ));
    const double error = error_vector->rms();

    energy_ = 0.5*(*aodensity_ * *hcore_).trace() + geom_->nuclear_repulsion();
    for (int i = 0; i != this->nocc(); ++i) energy_ += eig_[i];

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

    shared_ptr<Matrix1e> diis_density;
    if (iter >= diis_start_) {
      shared_ptr<Matrix1e> tmp_fock = diis.extrapolate(make_pair(fock, error_vector));
      shared_ptr<Matrix1e> intermediate(new Matrix1e(*tildex_ % *tmp_fock * *tildex_));
      intermediate->diagonalize(eig());
      shared_ptr<Coeff> tmp_coeff(new Coeff(*tildex_**intermediate));
      diis_density = tmp_coeff->form_density_rhf(nocc_);
    } else {
      diis_density = new_density;
    }

    aodensity_ = diis_density;
  }
}


shared_ptr<Reference> UHF::conv_to_ref() const {
  assert(false);
  shared_ptr<Reference> out(new Reference(geom_, coeff(), energy(), hcore(), schwarz(), nocc(), 0, geom_->nbasis()-nocc()));
  vector<double> e(eig_.get(), eig_.get()+geom_->nbasis());
  out->set_eig(e);
  return out;
}

