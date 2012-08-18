//
// Newint - Parallel electron correlation program.
// Filename: rysint.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


// Base class for the Rys-type integral evaluator

#ifndef __src_rysint_rysint_h
#define __src_rysint_rysint_h

#include <src/rysint/vrrlist.h>
#include <src/rysint/hrrlist.h>
#include <src/rysint/sortlist.h>
#include <src/scf/shell.h>
#include <src/stackmem.h>
#include <tuple>
#include <vector>
#include <memory>

class RysInt {
  protected:
    // some basic info for integral evaluations
    bool swap01_, swap23_;
    bool swap0123_;
    double AB_[3], CD_[3];
    int amapping_[ANG_VRR_END * ANG_VRR_END * ANG_VRR_END];
    int cmapping_[ANG_VRR_END * ANG_VRR_END * ANG_VRR_END];
    double *p_, *q_;
    double *xp_, *xq_, *coeff_, *coeffy_;
    double *T_, *U_;
    unsigned int contsize_, primsize_;
    size_t size_block_, size_alloc_;
    int prim0size_, prim1size_, prim2size_, prim3size_;
    int cont0size_, cont1size_, cont2size_, cont3size_;
    int asize_, csize_, amax_, amin_, cmax_, cmin_, amax1_, cmax1_;
    double *buff_;
    double *bkup_;

    std::vector<std::shared_ptr<const Shell> > basisinfo_;
    bool spherical_;

    // information on how many derivatives you take
    // 0 for ERI, 1 for gradients, etc. Set to 0 in the constructor, and will be
    // over written in the constructor of a derived class 
    int deriv_rank_;
    int tenno_;

    std::shared_ptr<HRRList> hrr_; 
    std::shared_ptr<VRRListBase> vrr_; 
    std::shared_ptr<SortList> sort_;

    double *data_;
    double *data2_;
    unsigned int size_final_;

    /// info for Rys quadruture
    double *roots_; 
    double *weights_;
    int rank_;

    // for screening
    int* screening_;
    int screening_size_;

    // init functions
    void set_swap_info(const bool swap_bra_ket = false);
    void set_ab_cd();
    void set_prim_contsizes();
    std::tuple<int, int, int, int> set_angular_info();
    // virtual init functions. The default is for ERI, NAI and their derivatives.
    // should be overloaded in Slater-type integrals
    virtual void root_weight(const int ps) = 0;
    virtual void compute_ssss(const double thr) = 0;
    void allocate_arrays(const size_t);
    void allocate_data(const int, const int, const int, const int);
    size_t size_allocated_;

    // for deallocation
    double* stack_save_;
    double* stack_save2_;


    // contraction
    void perform_contraction_new_outer(const int, const double*, const int, const int, double*,
                 const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int,
                 const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int);
    void perform_contraction_new_inner(const int, const int, const double*, const int, const int, double*,
                 const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int,
                 const std::vector<std::vector<double> >&, const std::vector<int>&, const std::vector<int>&, const int);
    // contraction for 1-e integrals
    void perform_contraction(const int, const double*, const int, const int, double*,
                             const std::vector<std::vector<double> >&, const std::vector<std::pair<int, int> >&, const int,
                             const std::vector<std::vector<double> >&, const std::vector<std::pair<int, int> >&, const int);

    std::shared_ptr<StackMem> stack_;

  public:
    RysInt(const std::vector<std::shared_ptr<const Shell> >&);
    ~RysInt();

    virtual void compute() = 0;

    /// retrieve a batch of integrals
    const double* data() const { return data_; };
    const double* data2() const { return data2_; };
    bool data2_exists() const { return data2_ != NULL; };
    unsigned int data_size() const { return size_final_; };

    size_t size_block() const { return size_block_; };

    bool swap01() const { return swap01_; }; 
    bool swap23() const { return swap23_; }; 
    bool swap0123() const { return swap0123_; }; 
};

#endif

