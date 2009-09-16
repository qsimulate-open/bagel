//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_osint_overlapbatch_h
#define __src_osint_overlapbatch_h

#include <vector>
#include <src/osint/osint.h>
#include <src/rysint/hrrlist.h>
#include <boost/shared_ptr.hpp>

class OverlapBatch : public OSInt {
  protected:
    HRRList hrr_;
    void perform_VRR(double*);

  public: 
    OverlapBatch(const std::vector<boost::shared_ptr<Shell> >&);
    ~OverlapBatch();

    void compute();
};

#endif
