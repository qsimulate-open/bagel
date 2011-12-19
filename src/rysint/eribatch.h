//
// Author: Toru Shiozaki
// Date  : April 2009
//
#ifndef __src_rysint_eriprim_h
#define __src_rysint_eriprim_h 

#include <cassert>
#include <vector>
#include <src/rysint/int2d.h>
#include <src/rysint/rysint.h>
#include <src/rysint/macros.h>
#include <memory>

class ERIBatch : public RysInt {

  protected:
    bool swap01_, swap23_;
    bool no_transpose_;

    double *p_, *q_;
    double *xp_, *xq_, *coeff_;
    double *T_;
    double AB_[3], CD_[3];
    unsigned int contsize_, primsize_, size_alloc_; 
    int prim0size_, prim1size_, prim2size_, prim3size_;
    int cont0size_, cont1size_, cont2size_, cont3size_;

    int asize_, csize_, amax_, amin_, cmax_, cmin_, amax1_, cmax1_;
    int amapping_[ANG_VRR_END * ANG_VRR_END * ANG_VRR_END];
    int cmapping_[ANG_VRR_END * ANG_VRR_END * ANG_VRR_END];

    /// buffer and intermediate storage
    double *buff_;
    double *bkup_;

    void perform_VRR1();
    void perform_VRR2();
    void perform_VRR3();
    void perform_VRR4();
    void perform_VRR5();
    void perform_VRR6();
    void perform_VRR7();
    void perform_VRR8();
    void perform_VRR9();
    void perform_VRR10();
    void perform_VRR11();
    void perform_VRR12();
    void perform_VRR13();
    void perform_VRR();

    void perform_contraction_new_outer(const int, const double*, const int, const int, double*,
                 const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int,
                 const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int);
    void perform_contraction_new_inner(const int, const double*, const int, const int, double*,
                 const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int,
                 const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int);

    void perform_HRR(const int, const double*, const double*, const double*, double*);

    void sort_indices(double*, const double*, const int, const int, const int, const int, const int, const bool);

  public:
    
    // dummy will never used.
    ERIBatch(const std::vector<std::shared_ptr<Shell> >, const double max_density, const double dummy = 0.0, const bool dum = true);
    ~ERIBatch();

    /// compute a batch of integrals
    void compute(); 

};

#endif

