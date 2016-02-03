//
// BAGEL - Parallel electron correlation program.
// Filename: spinfreebase.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_SMITH_SPINFREEBASE_H
#define __SRC_SMITH_SPINFREEBASE_H

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <chrono>
#include <src/ci/fci/civec.h>
#include <src/util/vec.h>
#include <src/smith/denom.h>
#include <src/smith/multitensor.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {

template<typename DataType>
class SpinFreeMethod {
  protected:
    using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;

    IndexRange virt_;
    IndexRange active_;
    IndexRange closed_;
    IndexRange all_;
    IndexRange ci_;

    // TODO these are redundant.
    std::shared_ptr<const IndexRange> rvirt_;
    std::shared_ptr<const IndexRange> ractive_;
    std::shared_ptr<const IndexRange> rclosed_;
    std::shared_ptr<const IndexRange> rci_;

    std::shared_ptr<const SMITH_Info<DataType>> info_;

    std::shared_ptr<const MatType> coeff_;
    double e0_;
    std::vector<double> e0all_;
    double core_energy_;
    double energy_;

    std::shared_ptr<Tensor_<DataType>> v2_;
    std::shared_ptr<Tensor_<DataType>> f1_;
    std::shared_ptr<Tensor_<DataType>> h1_;

    // contains the current RDMs to be used in smith
    std::shared_ptr<Tensor_<DataType>> rdm0_;
    std::shared_ptr<Tensor_<DataType>> rdm1_;
    std::shared_ptr<Tensor_<DataType>> rdm2_;
    std::shared_ptr<Tensor_<DataType>> rdm3_;
    std::shared_ptr<Tensor_<DataType>> rdm4_;

    // contains all the RDMs (for multistate runs)
    std::shared_ptr<Vec<Tensor_<DataType>>> rdm0all_;
    std::shared_ptr<Vec<Tensor_<DataType>>> rdm1all_;
    std::shared_ptr<Vec<Tensor_<DataType>>> rdm2all_;
    std::shared_ptr<Vec<Tensor_<DataType>>> rdm3all_;
    std::shared_ptr<Vec<Tensor_<DataType>>> rdm4all_;
    // the function to set RDMs to rdm1_, rdm2_, etc
    void set_rdm(const int jst, const int ist);

    // rdm ci derivatives
    std::shared_ptr<Tensor_<DataType>> rdm0deriv_;
    std::shared_ptr<Tensor_<DataType>> rdm1deriv_;
    std::shared_ptr<Tensor_<DataType>> rdm2deriv_;
    std::shared_ptr<Tensor_<DataType>> rdm3deriv_;
    std::shared_ptr<Tensor_<DataType>> rdm4deriv_;

    std::shared_ptr<Tensor_<DataType>> sigma_;

    // the diagonal denominator
    std::vector<double> eig_;

    // init functions
    void feed_rdm_denom(std::shared_ptr<const MatType> fockact);
    void feed_rdm_deriv(std::shared_ptr<const MatType> fockact);

    // printing functions called from the solve function of a derived class
    static void print_iteration();
    static void print_iteration(const int i, const double en, const double err, const double tim, const int istate = -1);
    static void print_iteration(const bool noconv);

    // compute e0 which is defined as Trace(f(x,x), gamma(x,x))
    void compute_e0();

    // denominator objects
    std::shared_ptr<const Denom<DataType>> denom_;

    // update t from the residual and denominator (this function does not zero out).
    void update_amplitude(std::shared_ptr<Tensor_<DataType>> t, std::shared_ptr<const Tensor_<DataType>> r) const;
    void update_amplitude(std::shared_ptr<MultiTensor_<DataType>> t, std::shared_ptr<const MultiTensor_<DataType>> r) const;

    // utility function
    void loop_over(std::function<void(const Index&, const Index&, const Index&, const Index&)>) const;
    // initialize t2 and r amplitude
    std::shared_ptr<Tensor_<DataType>> init_amplitude() const;
    std::shared_ptr<Tensor_<DataType>> init_residual() const;

  public:
    SpinFreeMethod(std::shared_ptr<const SMITH_Info<DataType>> r);
    ~SpinFreeMethod();

    IndexRange& virt() { return virt_; }
    IndexRange& all() { return all_; }
    IndexRange& closed() { return closed_; }

    std::shared_ptr<const SMITH_Info<DataType>> info() const { return info_; }
    std::shared_ptr<const MatType> coeff() const { return coeff_; }

    double e0() const { return e0_; }
    double energy() const { return energy_; }

    virtual void solve() = 0;

    DataType dot_product_transpose(std::shared_ptr<const Tensor_<DataType>> r, std::shared_ptr<const Tensor_<DataType>> t2) const;
    DataType dot_product_transpose(std::shared_ptr<const MultiTensor_<DataType>> r, std::shared_ptr<const MultiTensor_<DataType>> t2) const;
};

template<> void SpinFreeMethod<double>::feed_rdm_denom(std::shared_ptr<const Matrix>);
template<> void SpinFreeMethod<double>::feed_rdm_deriv(std::shared_ptr<const Matrix>);
template<> void SpinFreeMethod<std::complex<double>>::feed_rdm_denom(std::shared_ptr<const ZMatrix>);
template<> void SpinFreeMethod<std::complex<double>>::feed_rdm_deriv(std::shared_ptr<const ZMatrix>);
extern template class SpinFreeMethod<double>;
extern template class SpinFreeMethod<std::complex<double>>;

}
}

#endif
#endif
