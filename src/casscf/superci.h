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
    void common_init();

    void grad_vc(const std::shared_ptr<Matrix1e> fock, std::shared_ptr<RotFile> sigma);
    void grad_va(const std::shared_ptr<QFile> h, const std::shared_ptr<QFile> qxr, std::shared_ptr<RotFile> sigma);
    void grad_ca(const std::shared_ptr<Matrix1e> fock, const std::shared_ptr<Matrix1e> fock_inact, std::shared_ptr<RDM<1> > rdm1,
                 std::shared_ptr<QFile> qxr, std::shared_ptr<RotFile> sigma);

    void compute_qxr(double* int1ext, std::shared_ptr<RDM<2> > rdm2, std::shared_ptr<QFile> qxr);

    std::shared_ptr<Coeff> update_coeff(const std::shared_ptr<Coeff>, std::vector<double>);

    void sigma_at_at_(const std::shared_ptr<RotFile> cc, std::shared_ptr<RotFile> sigma,
                      const std::shared_ptr<QFile> gaa, const std::shared_ptr<Matrix1e> f);

    std::shared_ptr<RotFile> const_denom(const std::shared_ptr<QFile> gaa, const std::shared_ptr<Matrix1e> f) const;

    void update_orbitals(std::shared_ptr<RotFile> rot);

  public:
    SuperCI(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom);
    SuperCI(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<SCF> ref);
    ~SuperCI();

    void compute();

};

#endif
