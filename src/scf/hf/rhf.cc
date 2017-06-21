//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rhf.cc
// Copyright (C) 2013 Toru Shiozaki
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
#include <src/scf/hf/rhf.h>
#include <src/scf/hf/fock.h>
#include <src/prop/multipole.h>
#include <src/prop/sphmultipole.h>
#include <src/scf/dhf/population_analysis.h>
#include <src/util/muffle.h>

using namespace bagel;
using namespace std;

BOOST_CLASS_EXPORT_IMPLEMENT(RHF)

RHF::RHF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
 : SCF_base(idata, geom, re, !idata->get<bool>("df",true)), dodf_(idata->get<bool>("df",true)), restarted_(false) {

  cout << indent << "*** RHF ***" << endl << endl;
  if (nocc_ != noccB_) throw runtime_error("Closed shell SCF was called with nact != 0");

  // TODO In the future, this should probably be moved to SCF_base and designed to work properly there
  lshift_ = idata->get<double>("levelshift", 0.0);
  if (lshift_ != 0.0) {
    cout << "  level shift : " << setprecision(3) << lshift_ << endl << endl;
    levelshift_ = make_shared<ShiftVirtual<DistMatrix>>(nocc_, lshift_);
  }
}


void RHF::compute() {
  Timer scftime;

  shared_ptr<const Matrix> previous_fock = hcore_;

  shared_ptr<const Matrix> aodensity_;

  shared_ptr<const DistMatrix> tildex = tildex_->distmatrix();
  shared_ptr<const DistMatrix> hcore = hcore_->distmatrix();
  shared_ptr<const DistMatrix> overlap = overlap_->distmatrix();
  shared_ptr<const DistMatrix> coeff;
  shared_ptr<const DistMatrix> aodensity;

  auto newgeom = make_shared<const Geometry>(*geom_, "yang");
  auto fmmK = make_shared<const FMM>(newgeom, 5/*ns*/, 10/*lmaxJ*/, 1.0e-11/*fmmthresh*/, 0.1/*ws*/, true/*exchange*/, 2/*lmaxK*/);
  if (!restarted_) {
    if (coeff_ == nullptr) {
      shared_ptr<const DistMatrix> fock = previous_fock->distmatrix();
      if (geom_->spherical()) {
        auto aden = make_shared<const AtomicDensities>(geom_);
        shared_ptr<const Matrix> focka;
        if (!dofmm_) {
          if (dodf_) {
            focka = make_shared<const Fock<1>>(geom_, hcore_, aden, schwarz_);
          } else {
            focka = make_shared<const Fock<0>>(geom_, hcore_, aden, schwarz_);
          }
        } else {
#if 0
          shared_ptr<const Matrix> tmpJ = fmm_->compute_Fock_FMM(aden);
          //shared_ptr<const Matrix> tmpK = fmm_->compute_K_ff_from_den(aden, overlap_);
          focka = make_shared<const Matrix>(*hcore_ + *tmpJ);
          //focka = make_shared<const Matrix>(*hcore_ + *tmpJ - *tmpK);
#else
          shared_ptr<const Matrix> tmpJ = fmm_->compute_Fock_FMM_J(aden);
          shared_ptr<const Matrix> tmpKnf = fmmK->compute_Fock_FMM_K(aden);
          //shared_ptr<const Matrix> tmpK = fmmK->compute_K_ff_from_den(aden, overlap_);
          focka = make_shared<const Matrix>(*hcore_ + *tmpJ + *tmpKnf);
          //focka = make_shared<const Matrix>(*hcore_ + *tmpJ + *tmpKnf - *tmpK);
#endif
        }
        fock = focka->distmatrix();
      }
      DistMatrix intermediate = *tildex % *fock * *tildex;
      intermediate.diagonalize(eig());
      coeff = make_shared<const DistMatrix>(*tildex * intermediate);
    } else {
      shared_ptr<const Matrix> focka;
      if (!dofmm_) {
        if (!dodf_) {
          aodensity_ = coeff_->form_density_rhf(nocc_);
          focka = make_shared<const Fock<0>>(geom_, hcore_, aodensity_, schwarz_);
        } else {
          focka = make_shared<const Fock<1>>(geom_, hcore_, nullptr, coeff_->slice(0, nocc_), do_grad_, true/*rhf*/);
        }
      } else {
        aodensity_ = coeff_->form_density_rhf(nocc_);
#if 0
        shared_ptr<const Matrix> tmpJ = fmm_->compute_Fock_FMM(aodensity_);
        shared_ptr<const Matrix> tmpK = fmm_->compute_K_ff(make_shared<const Matrix>(coeff_->slice(0, nocc_)), overlap_);
        focka = make_shared<const Matrix>(*hcore_ + *tmpJ - *tmpK);
#else
        shared_ptr<const Matrix> tmpJ = fmm_->compute_Fock_FMM_J(aodensity_);
        shared_ptr<const Matrix> tmpKnf = fmmK->compute_Fock_FMM_K(aodensity_);
        shared_ptr<const Matrix> tmpK = fmmK->compute_K_ff(make_shared<const Matrix>(coeff_->slice(0, nocc_)), overlap_);
        focka = make_shared<const Matrix>(*hcore_ + *tmpJ + *tmpKnf - *tmpK);
#endif
      }
      DistMatrix intermediate = *tildex % *focka->distmatrix() * *tildex;
      intermediate.diagonalize(eig());
      coeff = make_shared<const DistMatrix>(*tildex * intermediate);
    }
    coeff_ = make_shared<const Coeff>(*coeff->matrix());

    cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
    cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
    cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
    scftime.tick_print("SCF startup");
    cout << endl;
    cout << indent << "=== RHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;

    diis_ = make_shared<DIIS<DistMatrix>>(diis_size_);
  } else {
    coeff = coeff_->distmatrix();
  }

  aodensity_ = coeff_->form_density_rhf(nocc_);
  aodensity = aodensity_->distmatrix();

  // starting SCF iteration
  shared_ptr<const Matrix> densitychange = aodensity_;

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer pdebug(1);

#ifndef DISABLE_SERIALIZATION
    if (restart_) {
      stringstream ss; ss << "scf_" << iter;
      OArchive archive(ss.str());
      archive << static_cast<Method*>(this);
    }
#endif

    auto ocoeff = make_shared<const Matrix>(coeff_->slice(0, nocc_));
    auto sc = make_shared<const Matrix>(*overlap_ * *ocoeff);
    if (!dofmm_) {
      if (!dodf_) {
        previous_fock = make_shared<Fock<0>>(geom_, previous_fock, densitychange, schwarz_);
        mpi__->broadcast(const_pointer_cast<Matrix>(previous_fock)->data(), previous_fock->size(), 0);
      } else {
        previous_fock = make_shared<Fock<1>>(geom_, hcore_, nullptr, coeff_->slice(0, nocc_), do_grad_, true/*rhf*/);
      }
    } else {
#if 0
      shared_ptr<const Matrix> tmpJ = fmm_->compute_Fock_FMM(aodensity_);
      shared_ptr<const Matrix> tmpK = fmm_->compute_K_ff(make_shared<const Matrix>(coeff_->slice(0, nocc_)), overlap_);
      previous_fock = make_shared<const Matrix>(*hcore_ + *tmpJ - *tmpK);
      //shared_ptr<const Matrix> tmpJ = fmm_->compute_Fock_FMM(densitychange);
      //previous_fock = make_shared<const Matrix>(*previous_fock + *tmpJ - *tmpK);
#else
      shared_ptr<const Matrix> tmpJ = fmm_->compute_Fock_FMM_J(aodensity_);
      shared_ptr<const Matrix> tmpKnf = fmmK->compute_Fock_FMM_K(aodensity_);
      shared_ptr<const Matrix> tmpK = fmmK->compute_K_ff(make_shared<const Matrix>(coeff_->slice(0, nocc_)), overlap_);
      previous_fock = make_shared<const Matrix>(*hcore_ + *tmpJ + *tmpKnf - *tmpK);
#endif
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

    if (diis_ || iter >= diis_start_) {
      fock = diis_->extrapolate({fock, error_vector});
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
    } else {
      aodensity_ = coeff_->form_density_rhf(nocc_);
    }
    aodensity = aodensity_->distmatrix();
    pdebug.tick_print("Post process");
  }
  // by default we compute dipoles
  if (!geom_->external() && multipole_print_) {
    if (dodf_) aodensity_ = aodensity->matrix();
    Multipole mu(geom_, aodensity_, multipole_print_);
    mu.compute();
    if (dma_print_ > 0) {
      SphMultipole smu(geom_, aodensity_, dma_print_);
      smu.compute();
    }
  }

  // print out orbital populations, if needed
  if (idata_->get<bool>("pop", false)) {
    cout << "    * Printing out population analysis to rhf.log" << endl;
    Muffle muf ("rhf.log");
    population_analysis(geom_, *coeff_, overlap_);
  }
}


shared_ptr<const Reference> RHF::conv_to_ref() const {
  auto out = make_shared<Reference>(geom_, coeff(), nocc(), 0, coeff_->mdim()-nocc(), vector<double>{energy_});
  out->set_eig(eig_);
  return out;
}
