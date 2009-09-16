//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_osint_kineticbatch_h
#define __src_osint_kineticbatch_h

#include <vector>
#include <src/osint/osint.h>
#include <boost/shared_ptr.hpp>

class KineticBatch : public OSInt {
  protected:
    void perform_VRR(double*);

  public: 
    KineticBatch(const std::vector<boost::shared_ptr<Shell> >&);
    ~KineticBatch();

    void compute();
};

#endif
