//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#ifndef __NEWINT_CASSCF_SUPERCI_H
#define __NEWINT_CASSCF_SUPERCI_H

#include <src/scf/scf.h>
#include <src/casscf/casscf.h>

class SuperCI : public CASSCF {

  protected:
    void update_orbitals_();
    void update_civectors_();
    void common_init();

  public:
    SuperCI(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom);
    SuperCI(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<SCF> ref);
    ~SuperCI();

    void compute();

};

#endif
