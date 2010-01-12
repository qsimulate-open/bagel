//
// Author: Toru Shiozaki
// Date  : June 2009
//
#ifndef __src_slater_slaterbatch_h
#define __src_slater_slaterbatch_h 

#include <cassert>
#include <vector>
#include <src/rysint/int2d.h>
#include <src/rysint/rysint.h>
#include <src/rysint/macros.h>
#include <src/slater/svrrlist.h>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

class SlaterBatch : public RysInt {

  protected:
    SVRRList svrr_;

    bool yukawa_;

    bool swap01_, swap23_;

    double gamma_;

    double *p_, *q_;
    double *xp_, *xq_, *coeff_, *coeffy_;
    double *T_, *U_;
    double AB_[3], CD_[3];
    unsigned int contsize_, primsize_, size_alloc_; 
    int prim0size_, prim1size_, prim2size_, prim3size_;
    int cont0size_, cont1size_, cont2size_, cont3size_;

    int asize_, csize_, amax_, amin_, cmax_, cmin_, amax1_, cmax1_;
    int amapping_[ANG_VRR_END * ANG_VRR_END * ANG_VRR_END];
    int cmapping_[ANG_VRR_END * ANG_VRR_END * ANG_VRR_END];

    std::vector<boost::tuple<int, double, double> > indexpair23_;

    /// buffer and intermediate storage
    double *buff_;
    double *bkup_;
    double *bkup2_;

    void perform_SVRR1();
    void perform_SVRR2();
    void perform_SVRR3();
    void perform_SVRR4();
    void perform_SVRR5();
    void perform_SVRR6();
    void perform_SVRR7();
    void perform_SVRR8();
    void perform_SVRR9();
    void perform_SVRR10();
    void perform_SVRR11();
    void perform_SVRR12();
    void perform_SVRR13();

    void perform_USVRR1();
    void perform_USVRR2();
    void perform_USVRR3();
    void perform_USVRR4();
    void perform_USVRR5();
    void perform_USVRR6();
    void perform_USVRR7();
    void perform_USVRR8();
    void perform_USVRR9();
    void perform_USVRR10();
    void perform_USVRR11();
    void perform_USVRR12();
    void perform_USVRR13();

    void perform_SVRR();
    void perform_USVRR();

    void perform_contraction_batch(const int, const double*, const int, const int, double*,
                                   const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int,
                                   const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int);
    void perform_contraction_new_outer(const int, const double*, const int, const int, double*,
                                       const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int,
                                       const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int);
    void perform_contraction_new_inner(const int, const double*, const int, const int, double*,
                                       const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int,
                                       const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int);

    void perform_HRR(const int, const double*, const double*, const double*, double*);

    void sort_indices(double*, const double*, const int, const int, const int, const int, const int, const bool);

  public:
    
    SlaterBatch(const std::vector<boost::shared_ptr<Shell> >, const double, const double gamma, const bool doyukawa = false);
    ~SlaterBatch();

    /// compute a batch of integrals
    void compute();

};

#endif

