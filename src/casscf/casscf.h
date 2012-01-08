//
// Author : Toru Shiozaki
// Date   : Dec 2011
// This is a base class for various CASSCF solvers.
// The assumpation is made that the CI and orbital optimizations are done independently.
// This should be a good one in a large systems.
//

#ifndef __NEWINT_CASSCF_CASSCF_H
#define __NEWINT_CASSCF_CASSCF_H

#include <cassert>
#include <vector>
#include <memory>
#include <src/scf/geometry.h>
#include <src/wfn/reference.h>
#include <src/fci/fci.h>
#include <src/fci/rdm.h>
#include <src/casscf/rotfile.h>

class CASSCF {

  protected:
    // input
    std::multimap<std::string, std::string> idata_; 

    // some internal information
    int nocc_; // sum of nact_ + nclosed_
    int nclosed_;
    int nact_;
    int nvirt_;
    int nbasis_;
    int nstate_;
    int max_iter_;
    int max_micro_iter_;
    double thresh_;
    double thresh_micro_;

    std::vector<double> occup_;
    std::shared_ptr<Coeff> coeff_natorb_;

    std::shared_ptr<Reference> ref_;
    std::shared_ptr<FCI> fci_;
    const std::shared_ptr<Geometry> geom_;
    void print_header() const;
    void common_init();

    void mute_stdcout();
    void resume_stdcout();

    std::shared_ptr<Matrix1e> ao_rdm1(std::shared_ptr<RDM<1> > rdm1, const bool inactive_only = false) const;

    std::shared_ptr<Fock<1> > hcore_;
    void one_body_operators(std::shared_ptr<Matrix1e>&, std::shared_ptr<QFile>&, std::shared_ptr<QFile>&, std::shared_ptr<QFile>&,
                            std::shared_ptr<RotFile>&, const bool superci=true);

    std::shared_ptr<Coeff> update_coeff(const std::shared_ptr<Coeff>, std::vector<double>) const;
    std::vector<double> form_natural_orbs();

  public:
    CASSCF(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<Reference> ref);
    virtual ~CASSCF();

    virtual void compute() { assert(false); };

    std::shared_ptr<Reference> ref() { return ref_; };
    std::shared_ptr<Reference> conv_to_ref() const {
      std::shared_ptr<Reference> out(new Reference(geom_, ref_->coeff(), ref_->hcore(), ref_->schwarz())); 
      return out;
    };

};

#endif
