//
// Author: Toru Shiozaki
// Date  : April 2009
//
#ifndef __src_rysint_naiprim_h
#define __src_rysint_naiprim_h 

#include <cassert>
#include <vector>
#include <src/scf/geometry.h>
#include <src/rysint/int2d.h>
#include <src/rysint/rysint.h>
#include <boost/shared_ptr.hpp>

class NAIBatch : public RysInt {

  protected:
    boost::shared_ptr<Geometry> geom_;
    int natom_;

    bool swap01_;

    double *p_;
    double *xp_, *coeff_;
    double *T_;
    double AB_[3];
    unsigned int contsize_, primsize_, size_alloc_; 
    int prim0size_, prim1size_;
    int cont0size_, cont1size_;

    int asize_, amax_, amin_, amax1_, asize_final_, asize_intermediate_;
    int amapping_[ANG_VRR_END * ANG_VRR_END * ANG_VRR_END];

    /// buffer and intermediate storage
    double *buff_;
    double *bkup_;

    /// for periodic calculations
    const int L_;
    const double A_;

    void perform_contraction(const int, const double*, const int, const int, double*,
                             const std::vector<std::vector<double> >&, const std::vector<std::pair<int, int> >&, const int,
                             const std::vector<std::vector<double> >&, const std::vector<std::pair<int, int> >&, const int);

  public:
    
    NAIBatch(const std::vector<boost::shared_ptr<Shell> >, const boost::shared_ptr<Geometry>, const int L = 0, const double A = 0.0);
    ~NAIBatch();

    /// compute a batch of integrals
    void compute();

};

#endif

