//
// BAGEL - Parallel electron correlation program.
// Filename: ks.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/ks/ks.h>
#include <src/prop/dipole.h>
#include <src/ks/dftgrid.h>
#include <src/ks/xcfunc.h>
#include <src/util/diis.h>

using namespace std;
using namespace bagel;

void KS::compute() {

  Matrix intermediate = *tildex_ % *hcore_ * *tildex_;
  intermediate.diagonalize(eig());
  coeff_ = shared_ptr<Coeff>(new Coeff(*tildex_ * intermediate));
  shared_ptr<Matrix> aodensity_ = coeff_->form_density_rhf(nocc_);

  Timer preptime;
  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;

  // TODO control from the input deck
//shared_ptr<DFTGrid_base> grid(new BLGrid(40, 194, geom_));
  shared_ptr<DFTGrid_base> grid(new DefaultGrid(geom_));
//shared_ptr<DFTGrid_base> grid(new TALGrid(75, 302, geom_));
//shared_ptr<DFTGrid_base> grid(new TALGrid(100, 770, geom_));
  preptime.tick_print("DFT grid generation");

  cout << indent << "     - DIIS with orbital gradients will be used." << endl << endl;
  cout << indent << "=== KS iteration (" << name_ << " / " << geom_->basisfile() << ") ===" << endl << indent << endl;

  DIIS<Matrix> diis(diis_size_);

  shared_ptr<Matrix> fock;

  Timer scftime;
  for (int iter = 0; iter != max_iter_; ++iter) {

    // fock operator without DFT xc 
    fock = shared_ptr<Matrix>(new Fock<1>(geom_, hcore_, aodensity_, coeff_->slice(0, nocc_), true, scale_ex_));

    // add xc 
    shared_ptr<const Matrix> xc;
    double exc;
    tie(xc, exc) = grid->compute_xc(name_, coeff_->slice(0, nocc_));

    energy_ = 0.5*((*hcore_+ *fock) * *aodensity_).trace() + exc + geom_->nuclear_repulsion();

    *fock += *xc;

    shared_ptr<const Matrix> error_vector(new Matrix(*fock**aodensity_**overlap_ - *overlap_**aodensity_**fock));

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
      fock = diis.extrapolate(make_pair(fock, error_vector));

    {
      shared_ptr<Matrix> intermediate(new Matrix(*tildex_ % *fock * *tildex_));
      intermediate->diagonalize(eig());
      coeff_ = shared_ptr<const Coeff>(new Coeff(*tildex_**intermediate));
    }
    aodensity_ = coeff_->form_density_rhf(nocc_);

  }

  // by default we compute dipoles
  if (!geom_->external()) {
    Dipole mu(geom_, aodensity_);
    mu.compute();
  }
}


shared_ptr<Reference> KS::conv_to_ref() const {
  shared_ptr<Reference> out(new Reference(geom_, coeff(), nocc(), 0, geom_->nbasis()-nocc(), energy()));
  vector<double> e(eig_.get(), eig_.get()+geom_->nbasis());
  out->set_eig(e);
  return out;
}

