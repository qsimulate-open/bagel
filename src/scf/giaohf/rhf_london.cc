//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rhf_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/scf/giaohf/rhf_london.h>
#include <src/scf/atomicdensities.h>
#include <src/wfn/zreference.h>
#include <src/prop/multipole.h>

using namespace bagel;
using namespace std;

BOOST_CLASS_EXPORT_IMPLEMENT(RHF_London)


RHF_London::RHF_London(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
 : SCF_base_London(idata, geom, re, !idata->get<bool>("df",true)), dodf_(idata->get<bool>("df",true)), restarted_(false) {

  cout << indent << "*** RHF ***" << endl << endl;
  if (nocc_ != noccB_) throw runtime_error("Closed shell SCF was called with nact != 0");

  // For the moment, I can't be bothered to test the level shifting apparatus for UHF and ROHF cases.
  // In the future, this should probably be moved to SCF_base and designed to work properly there
  lshift_ = idata->get<double>("levelshift", 0.0);
  if (lshift_ != 0.0) {
    cout << "  level shift : " << setprecision(3) << lshift_ << endl << endl;
    levelshift_ = make_shared<ShiftVirtual<DistZMatrix>>(nocc_, lshift_);
  }
}


void RHF_London::compute() {
  Timer scftime;

  shared_ptr<const ZMatrix> previous_fock = hcore_;
  shared_ptr<const ZMatrix> aodensity_;

  shared_ptr<const DistZMatrix> tildex = tildex_->distmatrix();
  shared_ptr<const DistZMatrix> tildextc = tildex_->transpose_conjg()->distmatrix();
  shared_ptr<const DistZMatrix> hcore = hcore_->distmatrix();
  shared_ptr<const DistZMatrix> overlap = overlap_->distmatrix();
  shared_ptr<const DistZMatrix> coeff;
  shared_ptr<const DistZMatrix> aodensity;

  if (!restarted_) {
    if (coeff_ == nullptr) {
      shared_ptr<const DistZMatrix> fock = hcore;
      if (dodf_ && geom_->spherical()) {
        auto aden = make_shared<const AtomicDensities>(geom_);
        auto zaden = std::make_shared<ZMatrix>(*aden, 1.0);
        auto focka = make_shared<const Fock_London<1>>(geom_, hcore_, zaden, schwarz_);
        fock = focka->distmatrix();
      }

      DistZMatrix intermediate = *tildex % *fock * *tildex;
      intermediate.diagonalize(eig());
      coeff = make_shared<const DistZMatrix>(*tildex * intermediate);
    } else {
      shared_ptr<const ZMatrix> focka;
      if (!dodf_) {
        aodensity_ = coeff_->form_density_rhf(nocc_, 0, 2.0);
        focka = make_shared<const Fock_London<0>>(geom_, hcore_, aodensity_, schwarz_);
      } else {
        focka = make_shared<const Fock_London<1>>(geom_, hcore_, nullptr, coeff_->slice(0, nocc_), do_grad_, true/*rhf*/);
      }
      DistZMatrix intermediate = *tildex % *focka->distmatrix() * *tildex;
      intermediate.diagonalize(eig());
      coeff = make_shared<const DistZMatrix>(*tildex * intermediate);
    }
    coeff_ = make_shared<const ZCoeff>(*coeff->matrix());

    cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
    cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
    cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
    scftime.tick_print("SCF startup");
    cout << endl;
    cout << indent << "=== RHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

    diis_ = make_shared<DIIS<DistZMatrix,ZMatrix>>(diis_size_);
  } else {
    coeff = coeff_->distmatrix();
  }

  if (!dodf_) {
    aodensity_ = coeff_->form_density_rhf(nocc_, 0, 2.0);
    aodensity = aodensity_->distmatrix();
  } else {
    aodensity = coeff->form_density_rhf(nocc_, 0, 2.0);
  }

  // starting SCF iteration
  shared_ptr<const ZMatrix> densitychange = aodensity_;

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer pdebug(1);

#ifndef DISABLE_SERIALIZATION
    if (restart_) {
      stringstream ss; ss << "scf_" << iter;
      OArchive archive(ss.str());
      archive << static_cast<Method*>(this);
    }
#endif
    if (!dodf_) {
      shared_ptr<ZMatrix> prev = previous_fock->copy();
      previous_fock = make_shared<Fock_London<0>>(geom_, previous_fock, densitychange, schwarz_);
      mpi__->broadcast(const_pointer_cast<ZMatrix>(previous_fock)->data(), previous_fock->size(), 0);
    } else {
      previous_fock = make_shared<Fock_London<1>>(geom_, hcore_, nullptr, coeff_->slice(0, nocc_), do_grad_, true/*rhf*/);
    }
    shared_ptr<const DistZMatrix> fock = previous_fock->distmatrix();

    const complex<double> zenergy  = 0.5*aodensity->dot_product(*hcore+*fock) + geom_->nuclear_repulsion();
    energy_  = real(zenergy);
    if (fabs(zenergy.imag()) > 1.0e-12) {
      stringstream ss; ss << "imaginary part of energy is nonzero!! Perhaps Fock is not Hermite for some reasons " << setprecision(10) << zenergy.imag();
//    throw runtime_error(ss.str());
      cout << ss.str() << endl;
    }

    pdebug.tick_print("Fock build");

    auto error_vector = make_shared<const DistZMatrix>(*fock**aodensity**overlap - *overlap**aodensity**fock);
    const double error = error_vector->rms();

    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                      << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;

    if (error < thresh_scf_) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      if (do_grad_) half_ = dynamic_pointer_cast<const Fock_London<1>>(previous_fock)->half();
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      break;
    }

    if (diis_ || iter >= diis_start_) {
      fock = diis_->extrapolate({fock, error_vector});
      pdebug.tick_print("DIIS");
    }

    DistZMatrix intermediate(*coeff % *fock * *coeff);

    if (levelshift_)
      levelshift_->shift(intermediate);

    intermediate.diagonalize(eig());

    pdebug.tick_print("Diag");

    coeff = make_shared<const DistZMatrix>(*coeff * intermediate);
    coeff_ = make_shared<const ZCoeff>(*coeff->matrix());

    if (!dodf_) {
      shared_ptr<const ZMatrix> new_density = coeff_->form_density_rhf(nocc_, 0, 2.0);
      densitychange = make_shared<ZMatrix>(*new_density - *aodensity_);
      aodensity_ = new_density;
      aodensity = aodensity_->distmatrix();
    } else {
      aodensity = coeff->form_density_rhf(nocc_, 0, 2.0);
    }
    pdebug.tick_print("Post process");
  }


  //TODO To report dipole moments with a GIAO basis, we'd need to implement multipole integrals over GIAO
#if 0
  if (!geom_->external()) {
    if (dodf_) aodensity_ = aodensity->matrix();
    Multipole mu(geom_, aodensity_, multipole_print_);
    mu.compute();
  }
#endif

}


shared_ptr<const Reference> RHF_London::conv_to_ref() const {
  auto out = make_shared<ZReference>(geom_, coeff_, energy_, nocc_, 0, coeff_->mdim()-nocc_);
  out->set_eig(eig_);
  return out;
}
