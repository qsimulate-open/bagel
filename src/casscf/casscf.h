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
#include <src/scf/scf.h>
#include <src/fci/fci.h>

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
    double thresh_;

    const std::shared_ptr<SCF> ref_;
    std::shared_ptr<FCI> fci_;
    const std::shared_ptr<Geometry> geom_;
    virtual void update_orbitals_() { assert(false); };
    virtual void update_civectors_() { assert(false); };
    void print_header() const;
    void common_init();

    void mute_stdcout();
    void resume_stdcout();

  public:
    CASSCF(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom);
    CASSCF(const std::multimap<std::string, std::string> idat, const std::shared_ptr<Geometry> geom, std::shared_ptr<SCF> ref);
    virtual ~CASSCF();

    virtual void compute() { assert(false); };

};

#endif
