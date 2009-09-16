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
#include <boost/shared_ptr.hpp>

class SCF {
  protected:
    const boost::shared_ptr<Geometry> geom_;
    const boost::shared_ptr<Overlap> overlap_;
    const boost::shared_ptr<Hcore> hcore_;
    boost::shared_ptr<TildeX> tildex_;
    boost::shared_ptr<Matrix1e> aodensity_;
    boost::shared_ptr<Coeff> coeff_;

    std::vector<double> shwarz_;
    void init_shwarz();

    double* eig_;
    void print_eig();

  public:
    SCF(const boost::shared_ptr<Geometry>);
    ~SCF();

    void compute();

    const boost::shared_ptr<Matrix1e> aodensity() { return aodensity_; };
};

#endif
