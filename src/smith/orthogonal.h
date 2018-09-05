//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: orthogonal.h
// Copyright (C) 2018 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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


#ifndef __SRC_SMITH_ORTHOGONAL_H
#define __SRC_SMITH_ORTHOGONAL_H

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

// class of excitations: note that airs and arsi are merged to one
enum Basis_Type { residual, amplitude };
enum Excitations { arbs, arbi, airj, risj, airs, arst, rist, aibj, total };

class Orthogonal_Basis {
  protected:

    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;

    std::vector<double> eig_;
    std::vector<double> e0all_;
    std::shared_ptr<Matrix> fockact_;

    // for intermediates
    std::vector<IndexRange> interm_;

    // orbital numbers
    size_t nact_;
    size_t nclosed_;
    size_t nvirt_;
    size_t nocc_;
    size_t ncore_;
    size_t nclo_;
    size_t norb_;
    size_t nstates_;

    bool sssr_;
    bool imag_;
    double shift_;

    double e0_;

    // rdms
    std::shared_ptr<Vec<Tensor_<double>>> rdm0all_;
    std::shared_ptr<Vec<Tensor_<double>>> rdm1all_;
    std::shared_ptr<Vec<Tensor_<double>>> rdm2all_;
    std::shared_ptr<Vec<Tensor_<double>>> rdm3all_;
    std::shared_ptr<Vec<Tensor_<double>>> rdm4all_;

    std::vector<size_t> size_;
    Basis_Type basis_type_;
    std::vector<std::shared_ptr<Matrix>> shalf_;

    std::vector<std::shared_ptr<MultiTensor_<double>>> data_;
    std::vector<std::shared_ptr<MultiTensor_<double>>> denom_;
    std::shared_ptr<Tensor_<double>> init_data(const int iext);
    std::shared_ptr<MultiTensor_<double>> get_contravariant(const int istate) const;
    void set_size(std::shared_ptr<const Denom<double>> d);
    void set_denom(std::shared_ptr<const Denom<double>> d);

    // copy of the functions in SpinFreeMethod
    std::shared_ptr<Tensor_<double>> init_amplitude() const;
    void loop_over(std::function<void(const Index&, const Index&, const Index&, const Index&)>) const;

    // feed rdm
    std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>,std::shared_ptr<RDM<3>>,std::shared_ptr<RDM<4>>> feed_rdm(const int ist, const int jst) const;

  public:
    // Orthogonal basis: construct from scratch (using denom)
    Orthogonal_Basis(std::shared_ptr<const SMITH_Info<double>> i, const IndexRange c, const IndexRange a, const IndexRange v, std::vector<double> f, std::vector<double> e0,
                     std::shared_ptr<const Matrix> fact, std::shared_ptr<const Denom<double>> d, const bool residual,
                     std::shared_ptr<Vec<Tensor_<double>>> g0, std::shared_ptr<Vec<Tensor_<double>>> g1, std::shared_ptr<Vec<Tensor_<double>>> g2,
                     std::shared_ptr<Vec<Tensor_<double>>> g3, std::shared_ptr<Vec<Tensor_<double>>> g4);
    // Or copy from existing orthogonal basis
    Orthogonal_Basis(const Orthogonal_Basis& o, const bool clone = true, const bool residual = true);

    // transform to orthogonal
    void transform_to_orthogonal(std::shared_ptr<const MultiTensor_<double>> t, const int istate);
    // Transform to redundant (should be amplitude)
    std::shared_ptr<MultiTensor_<double>> transform_to_redundant(const int istate) const;
    // update amplitude using residual
    void update(std::shared_ptr<const Orthogonal_Basis> residual, const int istate);
    // print convergence using source and residual
    void print_convergence(std::shared_ptr<const Orthogonal_Basis> source, std::shared_ptr<const Orthogonal_Basis> residual,
        const int istate, const int iter, const double error, const double timing) const;
    // add shift to residual using amplitude
    void add_shift(std::shared_ptr<const Orthogonal_Basis> amplitude, const int istate);
    // compute density matrix due to the shift
    std::tuple<std::shared_ptr<Matrix>,std::shared_ptr<Vec<double>>,
               std::shared_ptr<VecRDM<1>>,std::shared_ptr<VecRDM<2>>,std::shared_ptr<VecRDM<3>>,std::shared_ptr<VecRDM<3>>,std::vector<double>>
      make_d2_imag(std::shared_ptr<const Orthogonal_Basis> lambda) const;

    std::shared_ptr<MultiTensor_<double>>& data(const size_t i) { return data_[i]; }
    std::shared_ptr<MultiTensor_<double>> data(const size_t i) const { return data_[i]; }
    std::shared_ptr<MultiTensor_<double>> denom(const size_t i) const { return denom_[i]; }
    std::shared_ptr<Matrix> shalf(const int type) const { return shalf_[type]; }
    size_t size(const int type) const { return size_[type]; }
    size_t size_total() const { return size_[Excitations::total]; }
    std::vector<std::shared_ptr<VectorB>> vectorb() const;

    bool is_residual() const { return (basis_type_ == Basis_Type::residual); }
    bool is_amplitude() const { return (basis_type_ == Basis_Type::amplitude); }

    void zero(const int istate) { data_[istate]->zero(); }
    void zero() {
      for (auto& i : data_)
        i->zero();
    }
};

}
}

#endif
#endif
