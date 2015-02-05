
#ifndef __BAGEL_ASD_ORBOPT_BFGS_H
#define __BAGEL_ASD_ORBOPT_BFGS_H

#include <src/asd/orbopt/ras/rasscf.h>

namespace bagel {

class ASDRASBFGS : public ASDRASSCF {

  protected:
    std::shared_ptr<VectorB> preg_;

    void common_init() {
      std::cout << "    * Using the Quasi 2nd-order algorithm as noted in Chaban et al. TCA (1997)" << std::endl << std::endl;
    }

    // compute orbital gradients
    void grad_vc(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<ASDRASRotFile> sigma) const;
    void grad_va(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> qxr,   std::shared_ptr<Matrix> rdm1, std::shared_ptr<ASDRASRotFile> sigma) const;
    void grad_ca(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr, std::shared_ptr<Matrix> rdm1, std::shared_ptr<ASDRASRotFile> sigma) const;

    void grad_aa(std::shared_ptr<const Matrix> mcfock, std::shared_ptr<ASDRASRotFile> sigma) const;
    void grad_aa_with_preg(std::shared_ptr<const Matrix> mcfock, std::shared_ptr<ASDRASRotFile> sigma) const;

    // compute diagonal denominators
    std::shared_ptr<const ASDRASRotFile> compute_denom(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr, std::shared_ptr<const Matrix> rdm1, std::shared_ptr<const Matrix> mcfock) const;

  public:
    ASDRASBFGS(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr)
      : ASDRASSCF(idat, geom, ref) { 
      std::cout << "ASDRASBFGS(active-active) constructor" << std::endl; 
      common_init(); 
    }

    void compute() override;

};

}

#endif
