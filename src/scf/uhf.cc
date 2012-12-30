//
// BAGEL - Parallel electron correlation program.
// Filename: uhf.cc
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
#include <src/scf/uhf.h>
#include <src/prop/dipole.h>

using namespace std;
using namespace std::chrono;
using namespace bagel;

void UHF::compute() {

  string indent = "  ";

  if (coeff_ == nullptr || coeffB_ == nullptr) {
    Matrix intermediate = *tildex_ % *hcore_ * *tildex_;
    intermediate.diagonalize(eig());
    coeff_ = shared_ptr<Coeff>(new Coeff(*tildex_ * intermediate));
    coeffB_ = shared_ptr<Coeff>(new Coeff(*coeff_)); // since this is obtained with hcore
  }
  tie(aodensity_, aodensityA_, aodensityB_) = form_density_uhf();

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl;
  cout << endl;
  cout << indent << "    * DIIS with orbital gradients will be used."
            << endl << endl;
  cout << indent << "=== UHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  // starting SCF iteration
  eigB_ = unique_ptr<double[]>(new double[coeff_->mdim()]);

  DIIS<Matrix> diis(5);
  DIIS<Matrix> diisB(5);
  Timer scftime;
  for (int iter = 0; iter != max_iter_; ++iter) {

    shared_ptr<const Matrix> fockA(new Fock<1>(geom_, hcore_, aodensity_, schwarz_, coeff_->slice(0, nocc_)));
    shared_ptr<const Matrix> fockB(new Fock<1>(geom_, hcore_, aodensity_, schwarz_, coeffB_->slice(0, noccB_)));

    energy_ = 0.5*(*aodensity_ * *hcore_+ (*fockA * *aodensityA_ + *fockB * *aodensityB_)*0.5).trace() + geom_->nuclear_repulsion();

    shared_ptr<const Matrix> error_vector(new Matrix(*fockA**aodensityA_**overlap_ - *overlap_**aodensityA_**fockA
                                                    +*fockB**aodensityB_**overlap_ - *overlap_**aodensityB_**fockB));

    const double error = error_vector->rms();

    auto tp2 = high_resolution_clock::now();
    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2)
                                      << scftime.tick() << endl;

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) {
      fockA = diis.extrapolate(make_pair(fockA, error_vector));
      fockB = diisB.extrapolate(make_pair(fockB, error_vector));
    }

    {
      shared_ptr<Matrix> intermediate(new Matrix(*tildex_ % *fockA * *tildex_));
      intermediate->diagonalize(eig());
      coeff_ = shared_ptr<const Coeff>(new Coeff(*tildex_**intermediate));
    }{
      shared_ptr<Matrix> intermediate(new Matrix(*tildex_ % *fockB * *tildex_));
      intermediate->diagonalize(eigB());
      coeffB_ = shared_ptr<const Coeff>(new Coeff(*tildex_**intermediate));
    }
    tie(aodensity_, aodensityA_, aodensityB_) = form_density_uhf();

  }
  // now computes S^2 to see the spin contamination
  print_S2("UHF");

  // by default we compute dipoles
  if (!geom_->external()) {
    Dipole mu(geom_, aodensity_);
    mu.compute();
  }
}


void UHF::print_S2(const string tag) const {
  const double S2exact = (nocc_-noccB_)*(nocc_-noccB_+2)*0.25;
  const double contam = noccB_ - (*aodensityA_**overlap_**aodensityB_**overlap_).trace() * 0.25;
  cout << "    * S^2 (" << tag << ") is " << setprecision(4) << fixed << S2exact + contam << endl << endl;
}


tuple<shared_ptr<Coeff>, int, vector<shared_ptr<RDM<1> > > > UHF::natural_orbitals() const {
  shared_ptr<Matrix> cinv(new Matrix(*coeff_));
  cinv->inverse(); // TODO unnecessary
  shared_ptr<Matrix> intermediate(new Matrix(*cinv * *aodensity_ ^ *cinv));
  *intermediate *= -1.0;
  unique_ptr<double[]> occup(new double[geom_->nbasis()]);
  intermediate->diagonalize(occup.get());

  shared_ptr<Matrix> amat(new Matrix(*intermediate % (*cinv * *aodensityA_ ^ *cinv) * *intermediate));
  shared_ptr<Matrix> bmat(new Matrix(*intermediate % (*cinv * *aodensityB_ ^ *cinv) * *intermediate));

  int nocc = 0;
  // TODO adjust?
  const double tiny = 1.0e-10;
  for (int i = 0; i != geom_->nbasis(); ++i)
    if (occup[i] < -tiny) ++nocc;

  shared_ptr<RDM<1> > r(new RDM<1>(nocc));
  shared_ptr<RDM<1> > ra(new RDM<1>(nocc));
  shared_ptr<RDM<1> > rb(new RDM<1>(nocc));
  r->zero();
  for (int i = 0; i != nocc; ++i) r->element(i,i) = occup[i] * (-1.0);

  for (int i = 0; i != nocc; ++i) {
    for (int j = 0; j != nocc; ++j) {
      ra->element(j,i) = amat->element(j,i) * 0.5;
      rb->element(j,i) = bmat->element(j,i) * 0.5;
    }
  }

  vector<shared_ptr<RDM<1> > > rdm = {r, ra, rb};
  shared_ptr<Coeff> natorb(new Coeff(*coeff_ * *intermediate));

  return make_tuple(natorb, nocc, rdm);
}


shared_ptr<Reference> UHF::conv_to_ref() const {
  shared_ptr<Coeff> natorb;
  int nocc;
  vector<shared_ptr<RDM<1> > > rdm1;
  tie(natorb, nocc, rdm1) = natural_orbitals();
  shared_ptr<Reference> out(new Reference(geom_, natorb, 0, nocc, geom_->nbasis()-nocc, energy(), rdm1));

  // set alpha and beta coeffs
  out->set_coeff_AB(coeff_, coeffB_);

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
}


tuple<shared_ptr<Matrix>, shared_ptr<Matrix>, shared_ptr<Matrix> > UHF::form_density_uhf() const {
  shared_ptr<Matrix> outA = coeff_->form_density_rhf(nocc_);
  shared_ptr<Matrix> outB = coeffB_->form_density_rhf(noccB_);
  shared_ptr<Matrix> out(new Matrix(*outA+*outB));
  *out *= 0.5;
  return make_tuple(out, outA, outB);
}
