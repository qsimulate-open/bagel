//
// Author : Toru Shiozaki
// Date   : July 2009
//

#ifndef __src_pscf_pfock_h
#define __src_pscf_pfock_h

#include <src/pscf/pgeometry.h>
#include <src/util/pmatrix1e.h>
#include <src/pscf/phcore.h>
#include <src/util/pcompfile.h>
#include <src/rysint/eribatch.h>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <complex>

class PFock : public PMatrix1e {
  protected:
    const boost::shared_ptr<PFock> previous_;
    const boost::shared_ptr<PMatrix1e> density_;
    void pfock_two_electron_part();

    std::vector<double> schwarz_;
    int S2_;
    bool direct_;

    boost::shared_ptr<PCompFile<ERIBatch> > file_; 

  public:
    PFock(const boost::shared_ptr<PGeometry>, const boost::shared_ptr<PFock>, 
          const boost::shared_ptr<PMatrix1e>, const std::vector<double>&, const int, const bool dir);
    PFock(const boost::shared_ptr<PGeometry>, const boost::shared_ptr<PFock>, 
          const boost::shared_ptr<PMatrix1e>, const std::vector<double>&, const int, const bool dir, boost::shared_ptr<PCompFile<ERIBatch> >);
    PFock(const boost::shared_ptr<PGeometry>, const boost::shared_ptr<PHcore>);
    ~PFock();

    const bool direct() const { return direct_; };

};

#endif

