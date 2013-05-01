//
// BAGEL - Parallel electron correlation program.
// Filename: scf.h
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


#ifndef __BAGEL_SRC_SCF_SCF_H
#define __BAGEL_SRC_SCF_SCF_H

#include <src/scf/scf_base.h>
#include <src/scf/levelshift.h>
#include <src/util/diis.h>
#include <src/prop/dipole.h>
#include <src/wfn/reference.h>
#include <iostream>
#include <iomanip>
#include <src/parallel/mpi_interface.h>
#include <src/util/timer.h>
#include <src/scf/atomicdensities.h>

namespace bagel {

template<int DF>
class SCF : public SCF_base {
  protected:
    std::shared_ptr<LevelShift<DistMatrix>> levelshift_;

  public:
    SCF(const boost::property_tree::ptree& idata_, const std::shared_ptr<const Geometry> geom,
        const std::shared_ptr<const Reference> re = std::shared_ptr<const Reference>())
      : SCF_base(idata_, geom, re, DF==0) {

      std::cout << indent << "*** RHF ***" << std::endl << std::endl;
      if (nocc_ != noccB_) throw std::runtime_error("Closed shell SCF was called with nact != 0");

      // For the moment, I can't be bothered to test the level shifting apparatus for UHF and ROHF cases.
      // In the future, this should probably be moved to SCF_base and designed to work properly there
      double lshift = idata_.get<double>("levelshift", 0.0);
      if (lshift != 0.0) {
        std::cout << "  level shift : " << std::setprecision(3) << lshift << std::endl << std::endl;
        levelshift_ = std::shared_ptr<LevelShift<DistMatrix>>(new ShiftVirtual<DistMatrix>(nocc_, lshift));
      }

    }

    void compute() override {
      Timer scftime;

      std::shared_ptr<const Matrix> previous_fock = hcore_;
      std::shared_ptr<const Matrix> aodensity_;

      std::shared_ptr<const DistMatrix> tildex = tildex_->distmatrix();
      std::shared_ptr<const DistMatrix> hcore = hcore_->distmatrix();
      std::shared_ptr<const DistMatrix> overlap = overlap_->distmatrix();
      std::shared_ptr<const DistMatrix> coeff;
      std::shared_ptr<const DistMatrix> aodensity;

      if (coeff_ == nullptr) {
        std::shared_ptr<const DistMatrix> fock = hcore;
        if (DF == 1 && geom_->spherical()) {
          std::shared_ptr<const Matrix> aden(new AtomicDensities(geom_));
          std::shared_ptr<const Matrix> focka(new Fock<DF>(geom_, hcore_, aden, schwarz_));
          fock = focka->distmatrix();
        }
        DistMatrix intermediate = *tildex % *fock * *tildex;
        intermediate.diagonalize(eig());
        coeff = std::shared_ptr<const DistMatrix>(new DistMatrix(*tildex * intermediate));
      } else {
        std::shared_ptr<const Matrix> focka;
        if (DF == 0) {
          aodensity_ = coeff_->form_density_rhf(nocc_);
          focka = std::shared_ptr<const Matrix>(new Fock<DF>(geom_, hcore_, aodensity_, schwarz_));
        } else {
          focka = std::shared_ptr<const Matrix>(new Fock<DF>(geom_, hcore_, std::shared_ptr<const Matrix>(), coeff_->slice(0, nocc_), true));
        }
        DistMatrix intermediate = *tildex % *focka->distmatrix() * *tildex;
        intermediate.diagonalize(eig());
        coeff = std::shared_ptr<const DistMatrix>(new DistMatrix(*tildex * intermediate));
      }
      coeff_ = std::shared_ptr<const Coeff>(new Coeff(*coeff->matrix()));

      if (DF == 0) {
        aodensity_ = coeff_->form_density_rhf(nocc_);
        aodensity = aodensity_->distmatrix(); 
      } else {
        aodensity = coeff->form_density_rhf(nocc_);
      }

      std::cout << indent << "=== Nuclear Repulsion ===" << std::endl << indent << std::endl;
      std::cout << indent << std::fixed << std::setprecision(10) << std::setw(15) << geom_->nuclear_repulsion() << std::endl << std::endl;
      std::cout << indent << "    * DIIS with orbital gradients will be used." << std::endl << std::endl;
      scftime.tick_print("SCF startup");
      std::cout << std::endl;
      std::cout << indent << "=== RHF iteration (" + geom_->basisfile() + ") ===" << std::endl << indent << std::endl;

      // starting SCF iteration

      DIIS<DistMatrix> diis(diis_size_);
      std::shared_ptr<const Matrix> densitychange = aodensity_;

      for (int iter = 0; iter != max_iter_; ++iter) {
        Timer pdebug(1);

        if (DF == 0) {
          previous_fock = std::shared_ptr<Matrix>(new Fock<DF>(geom_, previous_fock, densitychange, schwarz_));
          mpi__->broadcast(previous_fock->data(), previous_fock->size(), 0);
        } else {
          previous_fock = std::shared_ptr<Matrix>(new Fock<DF>(geom_, hcore_, std::shared_ptr<const Matrix>(), coeff_->slice(0, nocc_), true));
        }
        std::shared_ptr<const DistMatrix> fock = previous_fock->distmatrix();

        energy_  = 0.5*aodensity->ddot(*hcore+*fock) + geom_->nuclear_repulsion();

        pdebug.tick_print("Fock build");

        std::shared_ptr<const DistMatrix> error_vector(new DistMatrix(*fock**aodensity**overlap - *overlap**aodensity**fock));
        const double error = error_vector->rms();

        std::cout << indent << std::setw(5) << iter << std::setw(20) << std::fixed << std::setprecision(8) << energy_ << "   "
                                          << std::setw(17) << error << std::setw(15) << std::setprecision(2) << scftime.tick() << std::endl;

        if (error < thresh_scf_) {
          std::cout << indent << std::endl << indent << "  * SCF iteration converged." << std::endl << std::endl;
          break;
        } else if (iter == max_iter_-1) {
          std::cout << indent << std::endl << indent << "  * Max iteration reached in SCF." << std::endl << std::endl;
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

        coeff = std::shared_ptr<const DistMatrix>(new DistMatrix(*coeff * intermediate));
        coeff_ = std::shared_ptr<const Coeff>(new Coeff(*coeff->matrix()));


        if (DF == 0) {
          std::shared_ptr<const Matrix> new_density = coeff_->form_density_rhf(nocc_);
          densitychange = std::shared_ptr<Matrix>(new Matrix(*new_density - *aodensity_));
          aodensity_ = new_density;
          aodensity = aodensity_->distmatrix();
        } else {
          aodensity = coeff->form_density_rhf(nocc_); 
        }
        pdebug.tick_print("Post process");
      }
      // by default we compute dipoles
      if (!geom_->external()) {
        if (DF != 0) aodensity_ = aodensity->matrix();
        Dipole mu(geom_, aodensity_);
        mu.compute();
      }
    }

    std::shared_ptr<Reference> conv_to_ref() const {
      std::shared_ptr<Reference> out(new Reference(geom_, coeff(), nocc(), 0, geom_->nbasis()-nocc(), energy()));
      std::vector<double> e(eig_.get(), eig_.get()+geom_->nbasis());
      out->set_eig(e);
      return out;
    }

};

}

#endif
