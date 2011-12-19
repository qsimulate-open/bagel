//
// Author: Toru Shiozaki
// Date  : May 2009
//
#ifndef __scf_scf_h
#define __scf_scf_h

#include <src/scf/geometry.h>
#include <src/scf/overlap.h>
#include <src/scf/hcore.h>
#include <src/scf/tildex.h>
#include <src/scf/fock.h>
#include <src/scf/coeff.h>
#include <memory>

class SCF {
  protected:
    const std::shared_ptr<Geometry> geom_;
    const std::shared_ptr<Overlap> overlap_;
    const std::shared_ptr<Hcore> hcore_;
    std::shared_ptr<TildeX> tildex_;
    std::shared_ptr<Matrix1e> aodensity_;
    std::shared_ptr<Coeff> coeff_;

    std::vector<double> shwarz_;
    void init_shwarz();

    double* eig_;
    void print_eig();

  public:
    SCF(const std::shared_ptr<Geometry>);
    ~SCF();

    void compute();

    const std::shared_ptr<Matrix1e> aodensity() { return aodensity_; };
    const std::shared_ptr<Coeff> coeff() { return coeff_; };
    const std::shared_ptr<Hcore> hcore() { return hcore_; };
};

#endif
