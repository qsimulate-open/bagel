//
// Author : Toru Shiozaki
// Date   : July 2009
//

#ifndef __src_pscf_pscf_disk_h
#define __src_pscf_pscf_disk_h

#include <boost/shared_ptr.hpp>
#include <src/pscf/pscf.h>
#include <src/util/pmatrix1e.h>
#include <src/util/pcompfile.h>

class PSCF_DISK : public PSCF {
  protected:
    void store_ERI();

  public:
    PSCF_DISK(const boost::shared_ptr<PGeometry>);
    PSCF_DISK(const boost::shared_ptr<PGeometry>, boost::shared_ptr<PMatrix1e>);
    ~PSCF_DISK();


};

#endif
