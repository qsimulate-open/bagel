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
  shared_ptr<const DistZMatrix> sohcore = sohcore_->distmatrix();
  shared_ptr<const DistZMatrix> sotildex = sotildex_->distmatrix();
  shared_ptr<const DistZMatrix> sooverlap = sooverlap_->distmatrix();
  shared_ptr<const DistZMatrix> socoeff;

  if (socoeff_ == nullptr) {
    shared_ptr<const DistZMatrix> sofock = sohcore;
    shared_ptr<DistZMatrix> intermediate = make_shared<DistZMatrix>(*sotildex % *sofock * *sotildex);
    intermediate->diagonalize(soeig());
    socoeff = make_shared<DistZMatrix>(*sotildex * *intermediate);
  } else {
    aodensity_ = socoeff->form_density_rhf(2*nocc_);
    shared_ptr<const ZMatrix> sofock = make_shared<const SOFock>(geom_, sohcore, aodensity_);
    shared_ptr<DistZMatrix> intermediate = make_shared<DistZMatrix>(*sotildex % *sofock * *sotildex);
    intermediate->diagonalize(soeig());
    socoeff = make_shared<DistZMatrix>(*sotildex * *intermediate);
  }
  aodensity_ = socoeff->form_density_rhf(2*nocc_);
  socoeff_ = make_shared<ZMatrix>(*socoeff->matrix());
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

  shared_ptr<const DistZMatrix> socoeff;
  for (int iter = 0; iter != max_iter_; ++iter) {
    shared_ptr<const ZMatrix> sofock = make_shared<const SOFock> (geom_, sohcore_, make_shared<ZMatrix>(socoeff_->slice(0, 2*nocc_)));
    const complex<double> energy = 0.5 * ((*sohcore_ + *sofock) * *aodensity_).trace() + geom_->nuclear_repulsion();
    assert(energy.imag() < 1e-8);
    energy_ = energy.real();
    auto error_vector = make_shared<const DistZMatrix>(*sofock * *aodensity_ * *sooverlap_ - *sooverlap_ * *aodensity_ * *sofock);
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
      const double twoe_energy = energy_ - geom_->nuclear_repulsion() - onee_energy;
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

    shared_ptr<DistZMatrix> intermediate = make_shared<DistZMatrix>(*sotildex_ % *sofock * *sotildex_);
    intermediate->diagonalize(soeig());
    socoeff = make_shared<const DistZMatrix>(*sotildex_ * *intermediate);
    socoeff_ = make_shared<const ZMatrix>(*socoeff->matrix());
    aodensity_ = socoeff->form_density_rhf(2*nocc_);
  }
}

shared_ptr<const ZMatrix> SOSCF::sotildex() {
  shared_ptr<const DistMatrix> tildex = tildex_->distmatrix();
  auto sotildex = make_shared<DistZMatrix>(2 * tildex->ndim(), 2 * tildex->mdim());
  sotildex->zero();
  sotildex->add_real_block(complex<double>(1.0, 0.0), 0, 0, tildex->ndim(), tildex->mdim(), *tildex);
  sotildex->add_real_block(complex<double>(1.0, 0.0), tildex->ndim(), tildex->mdim(), tildex->ndim(), tildex->mdim(), *tildex);
  auto out = make_shared<const ZMatrix>(*sotildex->matrix());
  return out;
}

shared_ptr<const ZMatrix> SOSCF::sooverlap() {
  shared_ptr<const DistMatrix> overlap = overlap_->distmatrix();
  auto sooverlap = make_shared<DistZMatrix>(2 * overlap->ndim(), 2 * overlap->mdim());
  sooverlap->zero();
  sooverlap->add_real_block(complex<double>(1.0, 0.0), 0, 0, overlap->ndim(), overlap->mdim(), *overlap);
  sooverlap->add_real_block(complex<double>(1.0, 0.0), overlap->ndim(), overlap->mdim(), overlap->ndim(), overlap->mdim(), *overlap);
  auto out = make_shared<const ZMatrix>(*sooverlap->matrix());
  return out;
}
