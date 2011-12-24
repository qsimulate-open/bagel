//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#ifndef __NEWINT_CASSCF_SUPERCI_H
#define __NEWINT_CASSCF_SUPERCI_H

#include <src/casscf/casscf.h>

class SuperCI : public CASSCF {

  protected:
    void update_orbitals_();
    void update_civectors_();

  public:
    SuperCI(const std::shared_ptr<Geometry> geom);
    ~SuperCI();

    void compute();

};

#endif
