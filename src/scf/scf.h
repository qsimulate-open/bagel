//
// Author: Toru Shiozaki
// Date  : May 2009
//
#ifndef __scf_scf_h
#define __scf_scf_h

#include <memory>
#include <string>
#include <map>
#include <src/scf/geometry.h>
#include <src/scf/overlap.h>
#include <src/scf/hcore.h>
#include <src/scf/tildex.h>
#include <src/scf/fock.h>
#include <src/scf/coeff.h>

class SCF {
  protected:
    std::multimap<std::string, std::string> idata_;
    const std::shared_ptr<Geometry> geom_;
    const std::shared_ptr<Overlap> overlap_;
    const std::shared_ptr<Hcore> hcore_;
    std::shared_ptr<TildeX> tildex_;
    std::shared_ptr<Matrix1e> aodensity_;
    std::shared_ptr<Coeff> coeff_;

    int max_iter_;
    double thresh_overlap_;
    double thresh_scf_;
    bool density_change_;

    std::vector<double> schwarz_;
    void init_schwarz();

    double* eig_;
    void print_eig();

  public:
    SCF(std::multimap<std::string, std::string>& idata_, const std::shared_ptr<Geometry>);
    ~SCF();

    void compute();

    const std::shared_ptr<Matrix1e> aodensity() { return aodensity_; };
    const std::shared_ptr<Coeff> coeff() { return coeff_; };
    const std::shared_ptr<Hcore> hcore() { return hcore_; };
    const std::vector<double>& schwarz() const { return schwarz_; };
};

#endif
