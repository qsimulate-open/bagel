//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rysintegral.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


// Base class for the Rys-type integral evaluator - now a class template

#ifndef __SRC_INTEGRAL_RYS_RYSINTEGRAL_H
#define __SRC_INTEGRAL_RYS_RYSINTEGRAL_H

#include <tuple>
#include <src/molecule/shell.h>
#include <src/util/constants.h>
#include <src/util/parallel/resources.h>
#include <src/integral/integral.h>

namespace bagel {

template <typename DataType, Int_t IntType = Int_t::Standard>
class RysIntegral : public Integral_base<DataType> {
  protected:
    // some basic info for integral evaluations
    bool swap01_, swap23_;
    bool swap0123_;
    std::array<double,3> AB_, CD_;
    int amapping_[ANG_VRR_END * ANG_VRR_END * ANG_VRR_END];
    int cmapping_[ANG_VRR_END * ANG_VRR_END * ANG_VRR_END];
    DataType *P_, *Q_;
    double *xp_, *xq_;
    DataType *coeff_, *coeffy_;
    DataType *T_, *U_;
    unsigned int contsize_, primsize_;
    size_t size_block_, size_alloc_;
    int prim0size_, prim1size_, prim2size_, prim3size_;
    int cont0size_, cont1size_, cont2size_, cont3size_;
    int asize_, csize_, amax_, amin_, cmax_, cmin_, amax1_, cmax1_;
    DataType *buff_;
    DataType *bkup_;

    std::array<std::shared_ptr<const Shell>,4> basisinfo_;
    bool spherical1_;
    bool spherical2_;

    // information on how many derivatives you take
    // 0 for ERI, 1 for gradients, etc. Set to 0 in the constructor, and will be
    // over written in the constructor of a derived class
    int deriv_rank_;
    int tenno_;
    int breit_;

    DataType *data_;
    double *data2_;
    unsigned int size_final_;

    /// info for Rys quadruture
    DataType *roots_;
    DataType *weights_;
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
    virtual void allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) = 0;

    void allocate_arrays(const size_t ps);

    size_t size_allocated_;

    // for deallocation
    DataType* stack_save_;
    double* stack_save2_;


    // contraction
    void perform_contraction_new_outer(const int nsize, const DataType* prim, const int pdim0, const int pdim1, DataType* cont,
                     const std::vector<std::vector<double>>& coeff0, const std::vector<int>& upper0, const std::vector<int>& lower0, const int cdim0,
                     const std::vector<std::vector<double>>& coeff1, const std::vector<int>& upper1, const std::vector<int>& lower1, const int cdim1);

    void perform_contraction_new_inner(const int nsize, const int ac, const DataType* prim, const int pdim0, const int pdim1, DataType* cont,
                     const std::vector<std::vector<double>>& coeff0, const std::vector<int>& upper0, const std::vector<int>& lower0, const int cdim0,
                     const std::vector<std::vector<double>>& coeff1, const std::vector<int>& upper1, const std::vector<int>& lower1, const int cdim1);

    // contraction for 1-e integrals
    void perform_contraction(const int asize, const DataType* prim, const int pdim0, const int pdim1, DataType* cont,
                               const std::vector<std::vector<double>>& coeff0, const std::vector<std::pair<int, int>>& ranges0, const int cdim0,
                               const std::vector<std::vector<double>>& coeff1, const std::vector<std::pair<int, int>>& ranges1, const int cdim1);

    bool allocated_here_;
    std::shared_ptr<StackMem> stack_;

  public:

    RysIntegral(const std::array<std::shared_ptr<const Shell>,4>& info, std::shared_ptr<StackMem> stack)
     : basisinfo_(info), spherical1_(info[0]->spherical()), spherical2_(info[2]->spherical()), deriv_rank_(0), tenno_(0), breit_(0) {
      assert(spherical1_ == info[1]->spherical());
      assert(spherical2_ == info[3]->spherical());

      if (stack == nullptr) {
        stack_ = resources__->get();
        allocated_here_ = true;
      } else {
        stack_ = stack;
        allocated_here_ = false;
      }
    }

    RysIntegral(const std::array<std::shared_ptr<const Shell>,2>& info, std::shared_ptr<StackMem> stack)
     : RysIntegral({{ info[0], info[1], std::make_shared<const Shell>(info[0]->spherical()), std::make_shared<const Shell>(info[0]->spherical()) }}, stack) { }


    ~RysIntegral() {
      // TODO this is a little inconsistent; stack should be allocated in the constructor of this class

      stack_->release(size_allocated_, buff_);
      if (tenno_) stack_->release(size_alloc_, stack_save2_);
      stack_->release(size_alloc_, stack_save_);

      if (allocated_here_)
        resources__->release(stack_);
    }

    virtual void compute() = 0;

    /// retrieve a batch of integrals
    virtual DataType* data(const int i) override { assert(i == 0); return data_; }
    const DataType* data() const { return data_; }
    const double* data2() const { return data2_; }
    bool data2_exists() const { return data2_ != nullptr; }
    size_t data_size() const { return size_final_; }

    size_t size_block() const { return size_block_; }

    bool swap01() const { return swap01_; }
    bool swap23() const { return swap23_; }
    bool swap0123() const { return swap0123_; }
};

using RysInt = RysIntegral<double>;

}

extern template class bagel::RysIntegral<double,bagel::Int_t::Standard>;
extern template class bagel::RysIntegral<std::complex<double>,bagel::Int_t::London>;

#endif

