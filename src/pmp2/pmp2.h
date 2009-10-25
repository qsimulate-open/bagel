//
// Author : Toru Shiozaki
// Date   : August 2009
//

#ifndef __src_pmp2_pmp2_h
#define __src_pmp2_pmp2_h

#include <vector>
#include <boost/shared_ptr.hpp>
#include <src/pscf/pcoeff.h>
#include <src/pscf/pgeometry.h>
#include <src/util/pfile.h>
#include <src/util/pcompfile.h>
#include <src/rysint/eribatch.h>

class PMP2 {
  protected:
    // Shared pointer for geometry.
    const boost::shared_ptr<PGeometry> geom_;

    // Coefficients for MO.
    const boost::shared_ptr<PCoeff> coeff_;

    // Coefficients for CABS (OBS / auxiliary part respectively)
    boost::shared_ptr<PMatrix1e> cabs_obs_;
    boost::shared_ptr<PMatrix1e> cabs_aux_;

    //
    const std::vector<double> eig_;

    boost::shared_ptr<PCompFile<ERIBatch> > ao_eri_;
    boost::shared_ptr<PMOFile<std::complex<double> > > eri_ii_pp_;
    int nbasis_;
    int nfrc_;
    int nocc_act_;
    int nocc_;
    int nvir_;
    int ncabs_;
    size_t noovv_;

    std::pair<boost::shared_ptr<PMatrix1e>, boost::shared_ptr<PMatrix1e> > generate_CABS();
    std::pair<boost::shared_ptr<PMatrix1e>, boost::shared_ptr<PMatrix1e> > generate_hJ();

  public:
    PMP2(const boost::shared_ptr<PGeometry>, const boost::shared_ptr<PCoeff>,
         const std::vector<double>, boost::shared_ptr<PCompFile<ERIBatch> >);
    ~PMP2();

    void compute();
    void compute_conv_mp2();

};

#endif

