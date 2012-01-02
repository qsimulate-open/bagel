//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __NEWINT_SRC_SCF_SCF_H
#define __NEWINT_SRC_SCF_SCF_H

#include <src/scf/scf_base.h>
#include <src/util/diis.h>
#include <iostream>
#include <iomanip>

template<int DF>
class SCF : public SCF_base {

  public:
    SCF(std::multimap<std::string, std::string>& idata_, const std::shared_ptr<Geometry> geom)
      : SCF_base(idata_, geom) { };

    ~SCF() {};

    void compute() {
      const bool highest_level = geom_->level() == 0;
      std::string indent = "  ";
      for (int i = 0; i != geom_->level(); ++i) indent += "|";
      if (!highest_level) indent += "  ";
    
      std::shared_ptr<Fock<DF> > previous_fock;
      {
        std::shared_ptr<Fock<DF> > fock(new Fock<DF>(geom_, hcore_));
        previous_fock = fock;
       
        Matrix1e intermediate = *tildex_ % *fock * *tildex_;
        intermediate.diagonalize(eig_);
        Coeff new_coeff(*tildex_ * intermediate);
        std::shared_ptr<Matrix1e> guess_density(new Matrix1e(new_coeff.form_density_rhf()));
    
        aodensity_ = guess_density;
      }
    
      if (highest_level) {
        std::cout << indent << "=== Nuclear Repulsion===" << std::endl << indent << std::endl;
        std::cout << indent << std::fixed << std::setprecision(10) << std::setw(15) << geom_->nuclear_repulsion() << std::endl;
        std::cout << std::endl; 
      }
      std::cout << indent << "=== RHF iteration (" + geom_->basisfile() + ")===" << std::endl << indent << std::endl;
    
      // starting SCF iteration
    
      DIIS<Matrix1e> diis(5);
      std::shared_ptr<Matrix1e> densitychange = aodensity_; // assumes hcore guess...
    
      for (int iter = 0; iter != max_iter_; ++iter) {
        int start = ::clock();
    
        std::shared_ptr<Fock<DF> > fock(new Fock<DF>(geom_, previous_fock, densitychange, schwarz_));
        previous_fock = fock;
    
        Matrix1e intermediate = *tildex_ % *fock * *tildex_;
        intermediate.diagonalize(eig_);
        std::shared_ptr<Coeff> new_coeff(new Coeff((*tildex_) * intermediate));
        coeff_ = new_coeff;
        std::shared_ptr<Matrix1e> new_density(new Matrix1e(coeff_->form_density_rhf()));
    
        std::shared_ptr<Matrix1e> error_vector(new Matrix1e(
          density_change_ ? (*new_density - *aodensity_) : (*fock * *aodensity_ * *overlap_ - *overlap_ * *aodensity_ * *fock)
        ));
        const double error = error_vector->rms();
    
        double energy = (*aodensity_ * *hcore_).trace() + geom_->nuclear_repulsion();
        for (int i = 0; i != geom_->nocc() / 2; ++i) energy += eig_[i];
    
        int end = ::clock();
        std::cout << indent << std::setw(5) << iter << std::setw(20) << std::fixed << std::setprecision(8) << energy << "   "
                                          << std::setw(17) << error << std::setw(15) << std::setprecision(2)
                                          << (end - start) / static_cast<double>(CLOCKS_PER_SEC) << std::endl; 
    
        if (error < thresh_scf_) {
          std::cout << indent << std::endl << indent << "  * SCF iteration converged." << std::endl << std::endl;
          break;
        } else if (iter == max_iter_-1) {
          std::cout << indent << std::endl << indent << "  * Max iteration reached in SCF." << std::endl << std::endl;
          break;
        }
    
        std::shared_ptr<Matrix1e> diis_density = diis.extrapolate(make_pair(new_density, error_vector));
    
        std::shared_ptr<Matrix1e> dtmp(new Matrix1e(*diis_density - *aodensity_));
        densitychange = dtmp; 
        aodensity_ = diis_density;
      }
    };
};

#endif
