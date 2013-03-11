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
  shared_ptr<DFTGrid_base> becke(new BLGrid(40, 194, geom_));
  preptime.tick_print("DFT grid generation");

  cout << indent << "     - DIIS with orbital gradients will be used." << endl << endl;
  cout << indent << "=== KS iteration (" << name_ << " / " << geom_->basisfile() << ") ===" << endl << indent << endl;

  DIIS<Matrix> diis(diis_size_);


  Timer scftime;
  for (int iter = 0; iter != max_iter_; ++iter) {

    // fock operator without DFT xc 
    shared_ptr<Matrix> fock;
    if (scale_ex_ != 0.0) {
      fock = shared_ptr<Matrix>(new Fock<1>(geom_, hcore_, aodensity_, coeff_->slice(0, nocc_), true, scale_ex_));
    } else {
      throw logic_error("pure DFT not yet implemented");
    }

    // add xc 
    *fock += *becke->compute_xcmat(name_, coeff_->slice(0, nocc_));

    energy_ = 0.5*(*aodensity_ * *hcore_+ *fock * *aodensity_).trace() + geom_->nuclear_repulsion();

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
#if 0
  shared_ptr<Coeff> natorb;
  int nocc;
  vector<shared_ptr<RDM<1>>> rdm1;
  tie(natorb, nocc, rdm1) = natural_orbitals();
  shared_ptr<Reference> out(new Reference(geom_, natorb, 0, nocc, geom_->nbasis()-nocc, energy(), rdm1));

  // set alpha and beta coeffs
  out->set_coeff_AB(coeff_, coeffB_);
  out->set_nocc(nocc_, noccB_);

  // compute an energy weighted 1RDM and store
  vector<double> ea(eig_.get(), eig_.get()+nocc_);
  vector<double> eb(eigB_.get(), eigB_.get()+nocc_);
  shared_ptr<Matrix> erdm = coeff_->form_weighted_density_rhf(nocc_, ea);
  *erdm += *coeffB_->form_weighted_density_rhf(noccB_, eb);
  *erdm *= 0.5;
  out->set_erdm1(erdm);

  // this is just dummy...
  vector<double> e(eig_.get(), eig_.get()+geom_->nbasis());
  out->set_eig(e);
  return out;
#else
  return shared_ptr<Reference>();
#endif
}

