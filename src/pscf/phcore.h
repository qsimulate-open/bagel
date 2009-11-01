//
// Author : Toru Shiozaki
// Date   : July 2009
//

#ifndef __src_pscf_phcore_h
#define __src_pscf_phcore_h

#include <src/util/pmatrix1e.h>
#include <src/pscf/pgeometry.h>
#include <boost/shared_ptr.hpp>

class PHcore : public PMatrix1e {
  protected:
    void computebatch(const std::vector<boost::shared_ptr<Shell> >&, const int, const int, const int, const int);

    // unfortunately, my implementation of CABS duplicates atoms... That means that
    // the nuclear attraction potential needs to be halved in Hcore that includes CABS.
    const bool cabs_;
    const bool kinetic_only_;

  public:
    PHcore(const boost::shared_ptr<PGeometry>, const bool cabs = false, const bool kinetic_only = false);
    ~PHcore();

};

#endif

