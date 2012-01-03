//
// Author: Toru Shiozaki
// Date  : May 2009
//
#ifndef __scf_scf_base_h
#define __scf_scf_base_h

#include <memory>
#include <string>
#include <map>
#include <src/scf/geometry.h>
#include <src/scf/overlap.h>
#include <src/scf/hcore.h>
#include <src/scf/tildex.h>
#include <src/scf/fock.h>
#include <src/scf/coeff.h>
#include <src/wfn/reference.h>

class SCF_base {
  protected:
    std::multimap<std::string, std::string> idata_;
    const std::shared_ptr<Geometry> geom_;
    const std::shared_ptr<Overlap> overlap_;
    const std::shared_ptr<Hcore> hcore_;
    std::shared_ptr<TildeX> tildex_;
    std::shared_ptr<Matrix1e> aodensity_;
    std::shared_ptr<Coeff> coeff_;

    int max_iter_;
    int diis_start_;
    double thresh_overlap_;
    double thresh_scf_;
    bool density_change_;

    std::vector<double> schwarz_;
    void init_schwarz();

    double* eig_;

  public:
    SCF_base(std::multimap<std::string, std::string>& idata_, const std::shared_ptr<Geometry>);
    ~SCF_base();

    virtual void compute() = 0;

    const std::shared_ptr<Geometry> geom() { return geom_; };
    const std::shared_ptr<Matrix1e> aodensity() { return aodensity_; };
    const std::shared_ptr<Coeff> coeff() { return coeff_; };
    void set_coeff(const std::shared_ptr<Coeff> o) { coeff_ = o; };
    const std::shared_ptr<Hcore> hcore() { return hcore_; };
    const std::vector<double>& schwarz() const { return schwarz_; };

    std::shared_ptr<Reference> conv_to_ref() {
      std::shared_ptr<Reference> out(new Reference(geom_, coeff(), hcore(), schwarz()));
      return out;
    };
};

#endif
