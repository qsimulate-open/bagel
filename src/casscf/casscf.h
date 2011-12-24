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
#include <src/scf/geometry.h>
#include <src/scf/scf.h>

class CASSCF {
  protected:
    const std::shared_ptr<SCF> ref_;
    const std::shared_ptr<Geometry> geom_;
    virtual void update_orbitals_() { assert(false); };
    virtual void update_civectors_() { assert(false); };
    void print_header() const;

  public:
    CASSCF(const std::shared_ptr<Geometry> geom);
    virtual ~CASSCF();

    virtual void compute() { assert(false); };

};

#endif
