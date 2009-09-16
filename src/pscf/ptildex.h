//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_pscf_ptildex_h
#define __src_pscf_ptildex_h

#include <complex>
#include <src/pscf/poverlap.h>
#include <src/util/pmatrix1e.h>
#include <boost/shared_ptr.hpp>

class PTildeX : public PMatrix1e {
  protected:

  public:
    PTildeX(const boost::shared_ptr<POverlap>);
    ~PTildeX();

};

#endif
