//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef _NEWINT_WFN_REFERENCE_H
#define _NEWINT_WFN_REFERENCE_H

#include <memory>
#include <src/scf/scf.h>
#include <src/scf/coeff.h>
#include <src/scf/geometry.h>

class Reference {
  protected:
    std::shared_ptr<Geometry> geom_; 
    std::shared_ptr<Coeff> coeff_;
    std::shared_ptr<Hcore> hcore_;
    std::vector<double> schwarz_;


  public:
    Reference(SCF& a);

    ~Reference() {};

    std::vector<double> schwarz() { return schwarz_; };
    std::shared_ptr<Hcore> hcore() { return hcore_; };
    const std::shared_ptr<Coeff> coeff() { return coeff_; };
    void set_coeff(const std::shared_ptr<Coeff> c) { coeff_ = c; };

};

#endif
