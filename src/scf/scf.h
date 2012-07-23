//
// Newint - Parallel electron correlation program.
// Filename: scf.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __NEWINT_SRC_SCF_SCF_H
#define __NEWINT_SRC_SCF_SCF_H

#include <src/scf/scf_base.h>
#include <src/util/diis.h>
#include <src/prop/dipole.h>
#include <iostream>
#include <iomanip>

template<int DF>
class SCF : public SCF_base {

  public:
    SCF(const std::multimap<std::string, std::string>& idata_, const std::shared_ptr<const Geometry> geom)
      : SCF_base(idata_, geom) {
      if (DF == 1) {
        // TODO init schwarz for auxiliary basis
      }
    };

    ~SCF() {};

    std::shared_ptr<Matrix1e> form_density_rhf() const { return coeff_->form_density_rhf(nocc_); };

    void compute() {
      const bool highest_level = geom_->level() == 0;
      std::string indent = "  ";
      for (int i = 0; i != geom_->level(); ++i) indent += "|";
      if (!highest_level) indent += "  ";
    
      std::shared_ptr<Fock<DF> > previous_fock;
      std::shared_ptr<Fock<DF> > hcore_fock;
      {
        previous_fock = std::shared_ptr<Fock<DF> >(new Fock<DF>(geom_, hcore_));
        if (DF) hcore_fock = previous_fock;
       
        Matrix1e intermediate = *tildex_ % *previous_fock * *tildex_;
        intermediate.diagonalize(eig());
        coeff_ = std::shared_ptr<Coeff>(new Coeff(*tildex_ * intermediate));
        aodensity_ = form_density_rhf();
      }
    
      if (highest_level) {
        std::cout << indent << "=== Nuclear Repulsion ===" << std::endl << indent << std::endl;
        std::cout << indent << std::fixed << std::setprecision(10) << std::setw(15) << geom_->nuclear_repulsion() << std::endl;
        std::cout << std::endl; 
      }
      std::cout << indent << "    * DIIS with " << (density_change_ ? "density changes" : "orbital gradients") << " will be used."
                << std::endl << std::endl;
      std::cout << indent << "=== RHF iteration (" + geom_->basisfile() + ") ===" << std::endl << indent << std::endl;
    
      // starting SCF iteration
    
      DIIS<Matrix1e> diis(5);
      std::shared_ptr<Matrix1e> densitychange = aodensity_; // assumes hcore guess...
    
      for (int iter = 0; iter != max_iter_; ++iter) {
        int start = ::clock();
    
        std::shared_ptr<Fock<DF> > fock;
        if (!DF) {
          fock = std::shared_ptr<Fock<DF> >(new Fock<DF>(geom_, previous_fock, densitychange, schwarz_));
        } else {
          fock = std::shared_ptr<Fock<DF> >(new Fock<DF>(geom_, hcore_fock, aodensity_, schwarz_));
        }
        previous_fock = fock;
    
        Matrix1e intermediate = *coeff_ % *fock * *coeff_;

// TODO level shift - needed?
//      intermediate.add_diag(1.0, this->nocc(), geom_->nbasis());

        intermediate.diagonalize(eig());
        coeff_ = std::shared_ptr<Coeff>(new Coeff((*coeff_) * intermediate));
        std::shared_ptr<Matrix1e> new_density = form_density_rhf();
    
        std::shared_ptr<Matrix1e> error_vector(new Matrix1e(
          density_change_ ? (*new_density - *aodensity_) : (*fock**aodensity_**overlap_ - *overlap_**aodensity_**fock)
        ));
        const double error = error_vector->rms();
    
        energy_ = 0.5*(*aodensity_ * *hcore_).trace() + geom_->nuclear_repulsion();
        for (int i = 0; i != this->nocc(); ++i) energy_ += eig_[i];
    
        int end = ::clock();
        std::cout << indent << std::setw(5) << iter << std::setw(20) << std::fixed << std::setprecision(8) << energy_ << "   "
                                          << std::setw(17) << error << std::setw(15) << std::setprecision(2)
                                          << (end - start) / static_cast<double>(CLOCKS_PER_SEC) << std::endl; 
    
        if (error < thresh_scf_) {
          std::cout << indent << std::endl << indent << "  * SCF iteration converged." << std::endl << std::endl;
          break;
        } else if (iter == max_iter_-1) {
          std::cout << indent << std::endl << indent << "  * Max iteration reached in SCF." << std::endl << std::endl;
          break;
        }
    
        std::shared_ptr<Matrix1e> diis_density;
        if (iter >= diis_start_) {
          std::shared_ptr<Matrix1e> tmp_fock = diis.extrapolate(make_pair(fock, error_vector));
          std::shared_ptr<Matrix1e> intermediate(new Matrix1e(*tildex_ % *tmp_fock * *tildex_));
          intermediate->diagonalize(eig());
          std::shared_ptr<Coeff> tmp_coeff(new Coeff(*tildex_**intermediate));
          diis_density = tmp_coeff->form_density_rhf(nocc_);
        } else {
          diis_density = new_density;
        }
    
        densitychange = std::shared_ptr<Matrix1e>(new Matrix1e(*diis_density - *aodensity_));
        aodensity_ = diis_density;
      }
      // by default we compute dipoles
      if (!geom_->external()) {
        Dipole mu(geom_, aodensity_);
        mu.compute();
      }
    };

    std::shared_ptr<Reference> conv_to_ref() const {
      std::shared_ptr<Reference> out(new Reference(geom_, coeff(), nocc(), 0, geom_->nbasis()-nocc(), energy()));
      std::vector<double> e(eig_.get(), eig_.get()+geom_->nbasis());
      out->set_eig(e);
      return out;
    };

};

#endif
