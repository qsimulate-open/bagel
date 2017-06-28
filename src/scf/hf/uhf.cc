//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: uhf.cc
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

#include <src/scf/atomicdensities.h>
#include <src/scf/hf/uhf.h>
#include <src/scf/hf/fock.h>
#include <src/prop/multipole.h>
#include <src/util/math/diis.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(UHF)

void UHF::initial_guess() {
  if (coeff_ == nullptr || coeffB_ == nullptr) {
    shared_ptr<const Matrix> fock = hcore_;
    if (geom_->spherical()) {
      auto aden = make_shared<const AtomicDensities>(geom_);
      fock = make_shared<const Fock<1>>(geom_, hcore_, aden, vector<double>());
    }
    Matrix intermediate = *tildex_ % *fock * *tildex_;
    intermediate.diagonalize(eig());
    coeff_ = make_shared<const Coeff>(*tildex_ * intermediate);
    coeffB_ = make_shared<const Coeff>(*coeff_);
  } else {
    tie(aodensity_, aodensityA_, aodensityB_) = form_density_uhf();
    auto fockA = make_shared<const Fock<1>>(geom_, hcore_, aodensity_, coeff_->slice(0, nocc_));
    auto fockB = make_shared<const Fock<1>>(geom_, hcore_, aodensity_, coeffB_->slice(0, noccB_));
    Matrix intermediateA = *tildex_ % *fockA * *tildex_;
    Matrix intermediateB = *tildex_ % *fockB * *tildex_;
    intermediateA.diagonalize(eig());
    intermediateB.diagonalize(eigB());
    coeff_ = make_shared<const Coeff>(*tildex_ * intermediateA);
    coeffB_ = make_shared<const Coeff>(*tildex_ * intermediateB);
  }
  tie(aodensity_, aodensityA_, aodensityB_) = form_density_uhf();
}

void UHF::compute() {

  initial_guess();

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl;
  cout << endl;
  cout << indent << "    * DIIS with orbital gradients will be used."
            << endl << endl;
  cout << indent << "=== UHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  // starting SCF iteration

  DIIS<Matrix> diis(diis_size_);
  DIIS<Matrix> diisB(diis_size_);
  Timer scftime;
  for (int iter = 0; iter != max_iter_; ++iter) {

    shared_ptr<const Matrix> fockA = make_shared<const Fock<1>>(geom_, hcore_, aodensity_, coeff_->slice(0, nocc_));
    shared_ptr<const Matrix> fockB = make_shared<const Fock<1>>(geom_, hcore_, aodensity_, coeffB_->slice(0, noccB_));

    energy_ = 0.25*((*hcore_+*fockA) * *aodensityA_ + (*hcore_+*fockB) * *aodensityB_).trace() + geom_->nuclear_repulsion();

    auto error_vector = make_shared<const Matrix>(*fockA**aodensityA_**overlap_ - *overlap_**aodensityA_**fockA
                                                 +*fockB**aodensityB_**overlap_ - *overlap_**aodensityB_**fockB);

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

    if (iter >= diis_start_) {
      fockA = diis.extrapolate({fockA, error_vector});
      fockB = diisB.extrapolate({fockB, error_vector});
    }

    {
      auto intermediate = make_shared<Matrix>(*tildex_ % *fockA * *tildex_);
      intermediate->diagonalize(eig());
      coeff_ = make_shared<const Coeff>(*tildex_**intermediate);
    }{
      auto intermediate = make_shared<Matrix>(*tildex_ % *fockB * *tildex_);
      intermediate->diagonalize(eigB());
      coeffB_ = make_shared<const Coeff>(*tildex_**intermediate);
    }
    tie(aodensity_, aodensityA_, aodensityB_) = form_density_uhf();

  }
  // now computes S^2 to see the spin contamination
  print_S2("UHF");

  // by default we compute dipoles
  if (!geom_->external() && multipole_print_) {
    Multipole mu(geom_, aodensity_, multipole_print_);
    scf_dipole_ = mu.compute();
  }
}


void UHF::print_S2(const string tag) const {
  const double S2exact = (nocc_-noccB_)*(nocc_-noccB_+2)*0.25;
  const double contam = noccB_ - (*aodensityA_**overlap_**aodensityB_**overlap_).trace() * 0.25;
  cout << "    * S^2 (" << tag << ") is " << setprecision(4) << fixed << S2exact + contam << endl << endl;
}


tuple<shared_ptr<Coeff>, int, shared_ptr<VecRDM<1>>> UHF::natural_orbitals() const {
  auto cinv = make_shared<Matrix>(*coeff_ % *overlap_);
  auto intermediate = make_shared<Matrix>(*cinv * *aodensity_ ^ *cinv);
  *intermediate *= -1.0;
  VectorB occup(coeff_->mdim());
  intermediate->diagonalize(occup);

  auto amat = make_shared<Matrix>(*intermediate % (*cinv * *aodensityA_ ^ *cinv) * *intermediate);
  auto bmat = make_shared<Matrix>(*intermediate % (*cinv * *aodensityB_ ^ *cinv) * *intermediate);

  int nocc = 0;
  // TODO adjust?
  const double tiny = 1.0e-10;
  for (int i = 0; i != coeff_->mdim(); ++i)
    if (occup[i] < -tiny) ++nocc;

  auto r = make_shared<RDM<1>>(nocc);
  auto ra = make_shared<RDM<1>>(nocc);
  auto rb = make_shared<RDM<1>>(nocc);
  r->zero();
  for (int i = 0; i != nocc; ++i) r->element(i,i) = occup[i] * (-1.0);

  for (int i = 0; i != nocc; ++i) {
    for (int j = 0; j != nocc; ++j) {
      ra->element(j,i) = amat->element(j,i) * 0.5;
      rb->element(j,i) = bmat->element(j,i) * 0.5;
    }
  }

  auto rdm = make_shared<VecRDM<1>>();
  rdm->emplace(0, r);
  rdm->emplace(1, ra);
  rdm->emplace(2, rb);

  return make_tuple(make_shared<Coeff>(*coeff_ * *intermediate), nocc, rdm);
}


shared_ptr<const Matrix> UHF::compute_erdm1() const {
  // compute an energy weighted 1RDM and store
  const VecView ea = eig_.slice(0, nocc_);
  const VecView eb = eigB_.slice(0, noccB_);
  shared_ptr<Matrix> erdm = coeff_->form_weighted_density_rhf(nocc_, ea);
  *erdm += *coeffB_->form_weighted_density_rhf(noccB_, eb);
  *erdm *= 0.5;
  return erdm;
}


shared_ptr<const Reference> UHF::conv_to_ref() const {
  shared_ptr<Coeff> natorb;
  int nocc;
  shared_ptr<VecRDM<1>> rdm1;
  tie(natorb, nocc, rdm1) = natural_orbitals();
  auto out = make_shared<Reference>(geom_, natorb, 0, nocc, coeff_->mdim()-nocc, vector<double>{energy_}, rdm1);

  // set alpha and beta coeffs
  out->set_coeff_AB(coeff_, coeffB_);
  out->set_nocc(nocc_, noccB_);

  // this is just dummy...
  out->set_eig(eig_);
  return out;
}


tuple<shared_ptr<const Matrix>, shared_ptr<const Matrix>, shared_ptr<const Matrix>> UHF::form_density_uhf() const {
  shared_ptr<const Matrix> outA = coeff_->form_density_rhf(nocc_);
  shared_ptr<const Matrix> outB = coeffB_->form_density_rhf(noccB_);
  auto out = make_shared<const Matrix>((*outA+*outB)*0.5);
  return make_tuple(out, outA, outB);
}
