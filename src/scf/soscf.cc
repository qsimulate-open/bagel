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
#include <src/rel/relreference.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(SOSCF)

SOSCF::SOSCF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
 : SCF_base(idata, geom, re), dodf_(idata->get<bool>("df",true)) {
  cout << indent << "*** Two-component ECP-SCF ***" << endl << endl;
  if (!dodf_)
    throw runtime_error("SOSCF requires density fitting!");

  soeig_ = VectorB(geom_->nbasis() * 2);
  sohcore_ = make_shared<SOHcore>(geom_, hcore_);
}

void SOSCF::compute() {
  Timer scftime;

  sooverlap_ = sooverlap();
  sotildex_ = sotildex();

  shared_ptr<const DistZMatrix> sohcore = sohcore_->distmatrix();
  shared_ptr<const DistZMatrix> sotildex = sotildex_->distmatrix();
  shared_ptr<const DistZMatrix> sooverlap = sooverlap_->distmatrix();
  shared_ptr<const DistZMatrix> socoeff, aodensity;

  if (socoeff_ == nullptr) {
    shared_ptr<const DistZMatrix> distfock = sohcore;
    auto intermediate = make_shared<DistZMatrix>(*sotildex % *distfock * *sotildex);
    intermediate->diagonalize(soeig());
    socoeff = make_shared<DistZMatrix>(*sotildex * *intermediate);
  } else {
    socoeff = socoeff_->distmatrix();
    aodensity = socoeff->form_density_rhf(2*nocc_);
    auto sofock = make_shared<const SOFock>(geom_, sohcore_, aodensity->matrix());
    shared_ptr<const DistZMatrix> distfock = sofock->distmatrix();
    auto intermediate = make_shared<DistZMatrix>(*sotildex % *distfock * *sotildex);
    intermediate->diagonalize(soeig());
    socoeff = make_shared<DistZMatrix>(*sotildex * *intermediate);
  }

  aodensity = socoeff->form_density_rhf(2*nocc_);

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
  scftime.tick_print("SOSCF startup");
  cout << endl;
  cout << indent << "=== SOSCF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

  DIIS<DistZMatrix, ZMatrix> diis(diis_size_);

  for (int iter = 0; iter != max_iter_; ++iter) {
    auto sofock = make_shared<const SOFock> (geom_, sohcore_, make_shared<ZMatrix>(socoeff->matrix()->slice(0, 2*nocc_)));
    shared_ptr<const DistZMatrix> distfock = sofock->distmatrix();
    const complex<double> energy = 0.5 * ((*sohcore_ + *sofock) * *aodensity->matrix()).trace() + geom_->nuclear_repulsion();
    assert(energy.imag() < 1e-8);
    energy_ = energy.real();
    auto error_vector = make_shared<const DistZMatrix>(*distfock * *aodensity * *sooverlap - *sooverlap * *aodensity * *distfock);
    const double error = error_vector->rms();

    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2) << scftime.tick();
    if (abs(energy.imag()) > 1e-12) {
      cout << "  *** Warning *** Im(E) = " << setw(15) << fixed << setprecision(12) << energy.imag() << endl;
    } else {
      cout << endl;
    }

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SOSCF iteration converged." << endl << endl;
      const double onee_energy = ((*sohcore_ * *aodensity->matrix()).trace()).real();
      const double twoe_energy = energy_ - geom_->nuclear_repulsion() - onee_energy;
      cout << indent << "    - One-electron energy" << setw(20) << fixed << setprecision(8) << onee_energy << endl;
      cout << indent << "    - Two-electron energy" << setw(20) << fixed << setprecision(8) << twoe_energy << endl;
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SOSCF." << endl << endl;
      break;
    }

    if (iter >= diis_start_) {
      distfock = diis.extrapolate({distfock, error_vector});
    }

    auto intermediate = make_shared<DistZMatrix>(*sotildex % *distfock * *sotildex);
    intermediate->diagonalize(soeig());
    auto socoeff = make_shared<const DistZMatrix>(*sotildex * *intermediate);
    aodensity = socoeff->form_density_rhf(2*nocc_);
  }

    socoeff_ = socoeff->matrix();
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

shared_ptr<const Reference> SOSCF::conv_to_ref() const {
  auto out =  make_shared<RelReference>(geom_, socoeff_, energy_, 0, 2*nocc_, 0, /*gaunt_*/ false, /*breit_*/ false);
  out->set_eig(soeig_);
  return out;
}
