//
// BAGEL - Parallel electron correlation program.
// Filename: scf.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/scf/scf.h>
#include <src/math/diis.h>
#include <src/prop/multipole.h>
#include <src/scf/atomicdensities.h>

using namespace bagel;
using namespace std;

SCF::SCF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
 : SCF_base(idata, geom, re, !idata->get<bool>("df",true)), dodf_(idata->get<bool>("df",true)) {

  cout << indent << "*** RHF ***" << endl << endl;
  if (nocc_ != noccB_) throw runtime_error("Closed shell SCF was called with nact != 0");

  // For the moment, I can't be bothered to test the level shifting apparatus for UHF and ROHF cases.
  // In the future, this should probably be moved to SCF_base and designed to work properly there
  double lshift = idata->get<double>("levelshift", 0.0);
  if (lshift != 0.0) {
    cout << "  level shift : " << setprecision(3) << lshift << endl << endl;
    levelshift_ = make_shared<ShiftVirtual<DistMatrix>>(nocc_, lshift);
  }
}


void SCF::compute() {
  Timer scftime;

  shared_ptr<const Matrix> previous_fock = hcore_;
  shared_ptr<const Matrix> aodensity_;

  shared_ptr<const DistMatrix> tildex = tildex_->distmatrix();
  shared_ptr<const DistMatrix> hcore = hcore_->distmatrix();
  shared_ptr<const DistMatrix> overlap = overlap_->distmatrix();
  shared_ptr<const DistMatrix> coeff;
  shared_ptr<const DistMatrix> aodensity;

  if (coeff_ == nullptr) {
    shared_ptr<const DistMatrix> fock = hcore;
    if (dodf_ && geom_->spherical()) {
      auto aden = make_shared<const AtomicDensities>(geom_);
      auto focka = make_shared<const Fock<1>>(geom_, hcore_, aden, schwarz_);
      fock = focka->distmatrix();
    }
    DistMatrix intermediate = *tildex % *fock * *tildex;
    intermediate.diagonalize(eig());
    coeff = make_shared<const DistMatrix>(*tildex * intermediate);
  } else {
    shared_ptr<const Matrix> focka;
    if (!dodf_) {
      aodensity_ = coeff_->form_density_rhf(nocc_);
      focka = make_shared<const Fock<0>>(geom_, hcore_, aodensity_, schwarz_);
    } else {
      focka = make_shared<const Fock<1>>(geom_, hcore_, shared_ptr<const Matrix>(), coeff_->slice(0, nocc_), do_grad_, true/*rhf*/);
    }
    DistMatrix intermediate = *tildex % *focka->distmatrix() * *tildex;
    intermediate.diagonalize(eig());
    coeff = make_shared<const DistMatrix>(*tildex * intermediate);
  }
  coeff_ = make_shared<const Coeff>(*coeff->matrix());

  if (!dodf_) {
    aodensity_ = coeff_->form_density_rhf(nocc_);
    aodensity = aodensity_->distmatrix();
  } else {
    aodensity = coeff->form_density_rhf(nocc_);
  }

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
  cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
  scftime.tick_print("SCF startup");
  cout << endl;
  cout << indent << "=== RHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  // starting SCF iteration

  DIIS<DistMatrix> diis(diis_size_);
  shared_ptr<const Matrix> densitychange = aodensity_;

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer pdebug(1);

    if (!dodf_) {
      previous_fock = make_shared<Fock<0>>(geom_, previous_fock, densitychange, schwarz_);
      mpi__->broadcast(const_pointer_cast<Matrix>(previous_fock)->data(), previous_fock->size(), 0);
    } else {
      previous_fock = make_shared<Fock<1>>(geom_, hcore_, shared_ptr<const Matrix>(), coeff_->slice(0, nocc_), do_grad_, true/*rhf*/);
    }
    shared_ptr<const DistMatrix> fock = previous_fock->distmatrix();

    energy_  = 0.5*aodensity->dot_product(*hcore+*fock) + geom_->nuclear_repulsion();

    pdebug.tick_print("Fock build");

    auto error_vector = make_shared<const DistMatrix>(*fock**aodensity**overlap - *overlap**aodensity**fock);
    const double error = error_vector->rms();

    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      if (do_grad_) half_ = dynamic_pointer_cast<const Fock<1>>(previous_fock)->half();
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) {
      fock = diis.extrapolate(make_pair(fock, error_vector));
      pdebug.tick_print("DIIS");
    }

    DistMatrix intermediate(*coeff % *fock * *coeff);

    if (levelshift_)
      levelshift_->shift(intermediate);

    intermediate.diagonalize(eig());
    pdebug.tick_print("Diag");

    coeff = make_shared<const DistMatrix>(*coeff * intermediate);
    coeff_ = make_shared<const Coeff>(*coeff->matrix());


    if (!dodf_) {
      shared_ptr<const Matrix> new_density = coeff_->form_density_rhf(nocc_);
      densitychange = make_shared<Matrix>(*new_density - *aodensity_);
      aodensity_ = new_density;
      aodensity = aodensity_->distmatrix();
    } else {
      aodensity = coeff->form_density_rhf(nocc_);
    }
    pdebug.tick_print("Post process");
  }
  // by default we compute dipoles
  if (!geom_->external()) {
    if (dodf_) aodensity_ = aodensity->matrix();
    Multipole mu(geom_, aodensity_, multipole_print_);
    mu.compute();
  }
}


shared_ptr<const Reference> SCF::conv_to_ref() const {
  auto out = make_shared<Reference>(geom_, coeff(), nocc(), 0, coeff_->mdim()-nocc(), energy());
  vector<double> e(eig_.get(), eig_.get()+coeff_->mdim());
  out->set_eig(e);
  return out;
}
