//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#ifndef __NEWINT_CASSCF_SUPERCI_H
#define __NEWINT_CASSCF_SUPERCI_H

#include <memory>
#include <string>
#include <map>
#include <src/scf/scf.h>
#include <src/casscf/casscf.h>
#include <src/casscf/rotfile.h>
#include <src/fci/rdm.h>

class SuperCI : public CASSCF {

  protected:
    void update_orbitals_();
    void update_civectors_();
    void common_init();

    void grad_vc(const std::shared_ptr<Matrix1e> fock, std::shared_ptr<RotFile> sigma);
    void grad_va(const std::shared_ptr<Matrix1e> fock_inact, std::shared_ptr<RDM<1> > rdm1,
                 std::shared_ptr<RotFile> qxr, std::shared_ptr<RotFile> sigma);
    void grad_ca(const std::shared_ptr<Matrix1e> fock, const std::shared_ptr<Matrix1e> fock_inact, std::shared_ptr<RDM<1> > rdm1,
                 std::shared_ptr<RotFile> qxr, std::shared_ptr<RotFile> sigma);

    void compute_qxr(double* int1ext, std::shared_ptr<RDM<2> > rdm2, std::shared_ptr<RotFile> qxr);

  public:
    SuperCI(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom);
    SuperCI(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<SCF> ref);
    ~SuperCI();

    void compute();

};

#endif
