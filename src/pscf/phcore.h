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

  public:
    PHcore(const boost::shared_ptr<PGeometry>);
    ~PHcore();

};

#endif

