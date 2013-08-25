//
// BAGEL - Parallel electron correlation program.
// Filename: dimer_scf.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Modified by: Shane Parker <shane.parker@u.northwestern.edu>
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


#include <src/dimer/dimer_scf.h>
#include <src/math/diis.h>
#include <src/prop/multipole.h>
#include <src/wfn/reference.h>

using namespace bagel;
using namespace std;

DimerSCF::DimerSCF(const std::shared_ptr<const PTree> idata, const shared_ptr<const Dimer> dimer)
  : SCF_base(idata, dimer->sgeom(), dimer->sref()) {

  dimer_ = dimer;

  const double shiftparameter = idata->get<double>("levelshift", 1.0e7);
  levelshift_ = make_shared<ShiftDimer>(dimer_, shiftparameter);
}

void DimerSCF::compute() {
  string indent = "  ";
  //aodensity_ = dimer_->form_density_rhf(coeff_);
  shared_ptr<Matrix> aodensity_ = coeff_->form_density_rhf(nocc_);

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
  cout << indent << "    * DIIS with orbital gradients will be used."
            << endl << endl;
  cout << indent << "=== DimerRHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  // starting SCF iteration

  DIIS<Matrix> diis(5);
  shared_ptr<Matrix> densitychange = aodensity_; // assumes hcore guess...

  //levelshift_->print_mo_data(coeff_);

  Timer timer;
  for (int iter = 0; iter != max_iter_; ++iter) {

    auto fock = make_shared<const Fock<1>>(geom_, hcore_, aodensity_, schwarz_);

    Matrix intermediate = *coeff_ % *fock * *coeff_;

    //intermediate.print("a", 24);

    levelshift_->shift(intermediate, coeff_);

    intermediate.diagonalize(eig());
    coeff_ = make_shared<Coeff>((*coeff_) * intermediate);
    shared_ptr<Matrix> new_density = dimer_->form_density_rhf(coeff_);

    auto error_vector = make_shared<const Matrix>(*fock**aodensity_**overlap_ - *overlap_**aodensity_**fock);
    const double error = error_vector->rms();

    // Note: the energy has to be computed this way since, by construction, the Fock matrix is NOT diagonal
    energy_ = 0.5*(*aodensity_ * (*fock + *hcore_)).trace() + geom_->nuclear_repulsion();
    //shared_ptr<Fock<1>> tmp_fock(new Fock<1>(geom_, hcore_fock, new_density, schwarz_));
    //energy_ = 0.5*(*new_density * (*tmp_fock + *hcore_)).trace() + geom_->nuclear_repulsion();

    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2) << timer.tick() << endl;

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      break;
    }

    // TODO - do I have anything to worry about while forming the diis_density? I think not but...
    shared_ptr<Matrix> diis_density;
    if (iter >= diis_start_) {
      shared_ptr<Matrix> tmp_fock = diis.extrapolate(make_pair(fock, error_vector));
      shared_ptr<Matrix> intermediate = make_shared<Matrix>(*tildex_ % *tmp_fock * *tildex_);
      intermediate->diagonalize(eig());
      shared_ptr<Coeff> tmp_coeff = make_shared<Coeff>(*tildex_**intermediate);
      diis_density = tmp_coeff->form_density_rhf(nocc_);
    } else {
      diis_density = new_density;
    }

    densitychange = make_shared<Matrix>(*diis_density - *aodensity_);
    aodensity_ = diis_density;
  }

  // by default we compute dipoles
  if (!geom_->external()) {
    Dipole mu(geom_, aodensity_);
    mu.compute();
  }

  //levelshift_->print_mo_data(coeff_);
}

shared_ptr<const Reference> DimerSCF::conv_to_ref() const {
  // Reorder here?
  auto out = make_shared<Reference>(geom_, coeff(), nocc(), 0, geom_->nbasis()-nocc(), energy());
  vector<double> e(eig_.get(), eig_.get()+geom_->nbasis());
  out->set_eig(e);
  return out;
}
