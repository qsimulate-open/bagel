//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef _NEWINT_WFN_REFERENCE_H
#define _NEWINT_WFN_REFERENCE_H

#include <memory>
#include <vector>
#include <src/scf/coeff.h>
#include <src/scf/hcore.h>
#include <src/scf/geometry.h>

class Reference {

  protected:
    std::shared_ptr<Geometry> geom_; 
    std::shared_ptr<Coeff> coeff_;
    std::shared_ptr<Hcore> hcore_;
    std::vector<double> schwarz_;

    const int nclosed_;
    const int nact_;
    const int nvirt_;

  public:
    Reference(std::shared_ptr<Geometry> g, std::shared_ptr<Coeff> c,
              std::shared_ptr<Hcore> h, const std::vector<double>& s,
              const int& nclo, const int& nact, const int& nvirt);

    ~Reference() {};

    std::shared_ptr<Geometry> geom() { return geom_; };
    std::vector<double> schwarz() { return schwarz_; };
    std::shared_ptr<Hcore> hcore() { return hcore_; };
    const std::shared_ptr<Coeff> coeff() { return coeff_; };
    void set_coeff(const std::shared_ptr<Coeff> c) { coeff_ = c; };

    int nclosed() const { return nclosed_; };
    int nact() const { return nact_; };
    int nvirt() const { return nvirt_; };

};

#endif
