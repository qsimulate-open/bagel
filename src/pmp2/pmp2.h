//
// Author : Toru Shiozaki
// Date   : August 2009
//

#ifndef __src_pmp2_pmp2_h
#define __src_pmp2_pmp2_h

#include <boost/shared_ptr.hpp>
#include <src/pscf/pcoeff.h>
#include <src/pscf/pgeometry.h>
#include <src/util/pfile.h>
#include <src/util/pcompfile.h>
#include <src/rysint/eribatch.h>

class PMP2 {
  protected:
    const boost::shared_ptr<PGeometry> geom_;
    const boost::shared_ptr<PCoeff> coeff_;
    const double* eig_;

    boost::shared_ptr<PCompFile<ERIBatch> > ao_eri_;
    boost::shared_ptr<PMOFile<std::complex<double> > > eri_ii_pp_;
    int nbasis_;
    int nfrc_;
    int nocc_act_;
    int nocc_;
    int nvir_;
    size_t noovv_;

  public:
    PMP2(const boost::shared_ptr<PGeometry>, const boost::shared_ptr<PCoeff>, const double*, boost::shared_ptr<PCompFile<ERIBatch> >);
    ~PMP2();

    void compute();
    void compute_conv_mp2();

};

#endif

