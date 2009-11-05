//
// Author : Toru Shiozaki
// Date   : August 2009
//

#ifndef __src_pmp2_pmp2_h
#define __src_pmp2_pmp2_h

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <src/pscf/pcoeff.h>
#include <src/pscf/pgeometry.h>
#include <src/util/pfile.h>
#include <src/util/pcompfile.h>
#include <src/util/pcompcabsfile.h>
#include <src/rysint/eribatch.h>
#include <src/slater/slaterbatch.h>

class PMP2 {
  typedef boost::shared_ptr<PMatrix1e> RefMatrix;
  typedef boost::shared_ptr<PCoeff> RefCoeff;
  typedef boost::shared_ptr<PMOFile<std::complex<double> > > RefMOFile;

  protected:
    // Shared pointer for geometry.
    const boost::shared_ptr<PGeometry> geom_;
    // PGeomery that has a union of OBS and CABS (in this order)
    // Will be initialized in generate_cabs()
    boost::shared_ptr<PGeometry> union_geom_;

    // Coefficients for MO.
    RefCoeff coeff_;

    // Coefficients for CABS (OBS / auxiliary part respectively)
    RefCoeff cabs_obs_;
    RefCoeff cabs_aux_;
    RefCoeff coeff_cabs_;
    RefMatrix coeff_entire_;

    // HF orbital energies
    const std::vector<double> eig_;

    // AO integrals used in several places
    boost::shared_ptr<PCompFile<ERIBatch> > eri_obs_;
    boost::shared_ptr<PCompFile<SlaterBatch> > stg_;
    boost::shared_ptr<PCompFile<SlaterBatch> > stg2_;
    boost::shared_ptr<PCompFile<SlaterBatch> > yp_;
    boost::shared_ptr<PCompCABSFile<ERIBatch> > eri_cabs_;
    boost::shared_ptr<PCompCABSFile<SlaterBatch> > stg_cabs_;
    boost::shared_ptr<PCompCABSFile<SlaterBatch> > stg_cabs2_;
    RefMOFile eri_ii_pp_;
    RefMOFile eri_ii_Ai_;
    RefMOFile stg_ii_pp_;
    RefMOFile stg_ii_Ai_;
    RefMOFile yp_ii_ii_;
    RefMOFile stg2_ii_ii_;

    // Target tensors
    RefMOFile X_;
    RefMOFile V_;
    RefMOFile B_;

    // Some misc constants
    int nbasis_;
    int nfrc_;
    int nocc_act_;
    int nocc_;
    int nvir_;
    int ncabs_;
    size_t noovv_;

    std::pair<RefCoeff, RefCoeff> generate_CABS();
    // Hartree and exchange matrix
    RefMatrix hJ_obs_obs_;
    RefMatrix hJ_obs_cabs_;
    RefMatrix hJ_cabs_obs_;
    RefMatrix hJ_cabs_cabs_;
    RefMatrix K_obs_obs_;
    RefMatrix K_obs_cabs_;
    RefMatrix K_cabs_obs_;
    RefMatrix K_cabs_cabs_;
    RefMatrix fock_obs_obs_;
    RefMatrix fock_obs_cabs_;
    RefMatrix fock_cabs_obs_;
    RefMatrix fock_cabs_cabs_;
    // Hartree builder
    const boost::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> generate_hJ() const;
    // Exchange builder
    const boost::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> generate_K() const;

    RefMatrix coulomb_runtime_OBS() const;
    RefMatrix exchange_runtime_OBS() const;
    RefMatrix coulomb_runtime() const;
    RefMatrix exchange_runtime() const;

  public:
    PMP2(const boost::shared_ptr<PGeometry>, const boost::shared_ptr<PCoeff>,
         const std::vector<double>, boost::shared_ptr<PCompFile<ERIBatch> >);
    ~PMP2();

    void compute();
    void compute_conv_mp2();

};

#endif

