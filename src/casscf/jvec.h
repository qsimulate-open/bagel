//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __NEWINT_SRC_CASSCF_JVEC_H
#define __NEWINT_SRC_CASSCF_JVEC_H

#include <src/fci/fci.h>
#include <src/scf/coeff.h>
#include <src/df/df.h>

class Jvec {
  protected:
    std::shared_ptr<DF_Half> half_; 
    std::unique_ptr<double[]> jvec_; 
    std::unique_ptr<double[]> rdm2_all_; 

  public:
    Jvec(std::shared_ptr<FCI> fci, std::shared_ptr<Coeff> c, const size_t, const size_t, const size_t);
    ~Jvec() {};

    const double* const half() const { return half_->data(); };
    const double* const jvec() const { return jvec_.get(); };
    const double* const rdm2_all() const { return rdm2_all_.get(); };

};

#endif
