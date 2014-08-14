//
// BAGEL - Parallel electron correlation program.
// Filename: soscf.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

#include <src/scf/soscf.h>
#include <src/scf/sofock.h>
#include <src/math/diis.h>
#include <src/scf/atomicdensities.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(SOSCF)

SOSCF::SOSCF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
 : SCF_base(idata, geom, re), dodf_(idata->get<bool>("df",true)) {
  cout << indent << "*** Two-component ECP-SCF ***" << endl << endl;
  if (!dodf_)
    throw runtime_error("SOSCF requires density fitting!");

  soeig_ = VectorB(geom_->nbasis() * 2);
  sohcore_base_ = make_shared<const SOHcore_base>(geom);
  sohcore_ = make_shared<SOHcore>(geom_, sohcore_base_);
}

void SOSCF::initial_guess() {
  sooverlap_ = sooverlap();
  sotildex_ = sotildex();

  if (socoeff_ == nullptr) {
    shared_ptr<const ZMatrix> sofock = sohcore_;
    shared_ptr<ZMatrix> intermediate = make_shared<ZMatrix>(*sotildex_ % *sofock * *sotildex_);
    intermediate->diagonalize(soeig());
    socoeff_ = make_shared<ZMatrix>(*sotildex_ * *intermediate);
  } else {
    aodensity_ = aodensity();
    auto sofock = make_shared<const SOFock>(geom_, sohcore_, make_shared<ZMatrix>(socoeff_->slice(0, nocc_ * 2)));
    shared_ptr<ZMatrix> intermediate = make_shared<ZMatrix>(*sotildex_ % *sofock * *sotildex_);
    intermediate->diagonalize(soeig());
    socoeff_ = make_shared<ZMatrix>(*sotildex_ * *intermediate);
  }
  aodensity_ = aodensity();
}

void SOSCF::compute() {
  Timer scftime;
  initial_guess();

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
  scftime.tick_print("SOSCF startup");
  cout << endl;
  cout << indent << "=== SOSCF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  DIIS<DistZMatrix, ZMatrix> diis(diis_size_);

  for (int iter = 0; iter != max_iter_; ++iter) {
    shared_ptr<const ZMatrix> sofock = make_shared<const SOFock> (geom_, sohcore_, make_shared<ZMatrix>(socoeff_->slice(0, nocc_ * 2)));
    const complex<double> energy = 0.5 * ((*sohcore_ + *sofock) * *aodensity_).trace() + geom_->nuclear_repulsion();
    assert(energy.imag() < 1e-8);
    energy_ = energy.real();
    auto error_vector = make_shared<const ZMatrix>(*sofock * *aodensity_ * *sooverlap_ - *sooverlap_ * *aodensity_ * *sofock);
    auto real_error_vector = error_vector->get_real_part();
    const double error = real_error_vector->rms();

    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2) << scftime.tick();
    if (abs(energy.imag()) > 1e-12) {
      cout << "  *** Warning *** Im(E) = " << setw(15) << fixed << setprecision(12) << energy.imag() << endl;
    } else {
      cout << endl;
    }

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SOSCF iteration converged." << endl << endl;
      const double onee_energy = ((*sohcore_ * *aodensity_).trace()).real();
      const double twoe_energy = energy_ - onee_energy;
      cout << indent << "    - One-electron energy" << setw(20) << fixed << setprecision(8) << onee_energy << endl;
      cout << indent << "    - Two-electron energy" << setw(20) << fixed << setprecision(8) << twoe_energy << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SOSCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) {
      sofock = diis.extrapolate({sofock, error_vector});
    }

    shared_ptr<ZMatrix> intermediate = make_shared<ZMatrix>(*sotildex_ % *sofock * *sotildex_);
    intermediate->diagonalize(soeig());
    socoeff_ = make_shared<ZMatrix>(*sotildex_ * *intermediate);
    aodensity_ = aodensity();
  }
}

shared_ptr<const ZMatrix> SOSCF::sotildex() {
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(2 * tildex_->ndim(), 2 * tildex_->mdim());
  out->zero();
  out->add_real_block(complex<double>(1.0, 0.0), 0, 0, tildex_->ndim(), tildex_->mdim(), *tildex_);
  out->add_real_block(complex<double>(1.0, 0.0), tildex_->ndim(), tildex_->mdim(), tildex_->ndim(), tildex_->mdim(), *tildex_);
  return out;
}

shared_ptr<const ZMatrix> SOSCF::sooverlap() {
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(2 * overlap_->ndim(), 2 * overlap_->mdim());
  out->zero();
  out->add_real_block(complex<double>(1.0, 0.0), 0, 0, overlap_->ndim(), overlap_->mdim(), *overlap_);
  out->add_real_block(complex<double>(1.0, 0.0), overlap_->ndim(), overlap_->mdim(), overlap_->ndim(), overlap_->mdim(), *overlap_);
  return out;
}

std::shared_ptr<const ZMatrix> SOSCF::aodensity() {
  const int n = geom_->nbasis();
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(2*n, 2*n);
  out->zero();

  shared_ptr<const ZMatrix> coeffa = socoeff_->get_submatrix(0, 0, n, nocc_*2);
  shared_ptr<const ZMatrix> coeffb = socoeff_->get_submatrix(n, 0, n, nocc_*2);

  shared_ptr<const Matrix> rcoeffa = coeffa->get_real_part();
  shared_ptr<const Matrix> icoeffa = coeffa->get_imag_part();
  shared_ptr<const Matrix> rcoeffb = coeffb->get_real_part();
  shared_ptr<const Matrix> icoeffb = coeffb->get_imag_part();

  shared_ptr<const Matrix> reDaa = make_shared<const Matrix>((*rcoeffa ^ *rcoeffa) + (*icoeffa ^ *icoeffa));
  shared_ptr<const Matrix> imDaa = make_shared<const Matrix>(((*rcoeffa ^ *icoeffa) * (-1.0)) + (*icoeffa ^ *rcoeffa));
  out->add_real_block(complex<double>(1.0, 0.0), 0, 0, n , n, *reDaa);
  out->add_real_block(complex<double>(0.0, 1.0), 0, 0, n , n, *imDaa);

  shared_ptr<const Matrix> reDbb = make_shared<const Matrix>((*rcoeffb ^ *rcoeffb) + (*icoeffb ^ *icoeffb));
  shared_ptr<const Matrix> imDbb = make_shared<const Matrix>(((*rcoeffb ^ *icoeffb) * (-1.0)) + (*icoeffb ^ *rcoeffb));
  out->add_real_block(complex<double>(1.0, 0.0), n, n, n , n, *reDbb);
  out->add_real_block(complex<double>(0.0, 1.0), n, n, n , n, *imDbb);

  shared_ptr<const Matrix> reDab = make_shared<const Matrix>((*rcoeffa ^ *rcoeffb) + (*icoeffa ^ *icoeffb));
  shared_ptr<const Matrix> imDab = make_shared<const Matrix>(((*rcoeffa ^ *icoeffb) * (-1.0)) + (*icoeffa ^ *rcoeffb));
  out->add_real_block(complex<double>(1.0, 0.0), 0, n, n , n, *reDab);
  out->add_real_block(complex<double>(0.0, 1.0), 0, n, n , n, *imDab);

  shared_ptr<const Matrix> reDba = make_shared<const Matrix>((*rcoeffb ^ *rcoeffa) + (*icoeffb ^ *icoeffa));
  shared_ptr<const Matrix> imDba = make_shared<const Matrix>(((*rcoeffb ^ *icoeffa) * (-1.0)) + (*icoeffb ^ *rcoeffa));
  out->add_real_block(complex<double>(1.0, 0.0), n, 0, n , n, *reDba);
  out->add_real_block(complex<double>(0.0, 1.0), n, 0, n , n, *imDba);

  return out;
}

