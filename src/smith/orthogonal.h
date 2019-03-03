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
// TODO need fixes
  protected:
    std::map<int,std::string> to_denom_;

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

    Basis_Type basis_type_;
    std::shared_ptr<const Denom<double>> d_;

    std::vector<std::shared_ptr<MultiTensor_<double>>> data_;
    std::vector<std::shared_ptr<MultiTensor_<double>>> denom_;
    std::shared_ptr<Tensor_<double>> init_data(const int iext, const int istate);
    std::shared_ptr<MultiTensor_<double>> weight_by_denom(const int istate, std::shared_ptr<const MultiTensor_<double>> original) const;
    void set_shalf(std::shared_ptr<const Denom<double>> d);
    void set_denom(std::shared_ptr<const Denom<double>> d);

    // copy of the functions in SpinFreeMethod
    std::shared_ptr<Tensor_<double>> init_amplitude() const;
    void loop_over(std::function<void(const Index&, const Index&, const Index&, const Index&)>) const;

  public:
    // Orthogonal basis: construct from scratch (using denom)
    Orthogonal_Basis(std::shared_ptr<const SMITH_Info<double>> i, const IndexRange c, const IndexRange a, const IndexRange v, std::vector<double> f, std::vector<double> e0,
                     std::shared_ptr<const Matrix> fact, std::shared_ptr<const Denom<double>> d, const bool residual);
    // Or copy from existing orthogonal basis
    Orthogonal_Basis(const Orthogonal_Basis& o, const bool clone = true, const bool residual = true);

    IndexRange& virt() { return virt_; }
    IndexRange& closed() { return closed_; }
    IndexRange& interm(const int dataindex) { return interm_[dataindex]; }

    // transform to orthogonal
    void transform_to_orthogonal(std::shared_ptr<const MultiTensor_<double>> t, const int istate);
    // Transform to redundant (should be amplitude)
    std::shared_ptr<MultiTensor_<double>> transform_to_redundant(const int istate) const;
    // update amplitude using residual
    void update(std::shared_ptr<const Orthogonal_Basis> residual, const int istate);
    // print convergence using source and residual
    double print_convergence(std::shared_ptr<const Orthogonal_Basis> source, std::shared_ptr<const Orthogonal_Basis> residual,
        const int istate, const int iter, const double error, const double timing) const;
    // add shift to residual using amplitude
    void add_shift(std::shared_ptr<const Orthogonal_Basis> amplitude, const int istate);

    std::shared_ptr<MultiTensor_<double>>& data(const size_t i) { return data_[i]; }
    std::shared_ptr<MultiTensor_<double>> data(const size_t i) const { return data_[i]; }
    std::shared_ptr<MultiTensor_<double>> denom(const size_t i) const { return denom_[i]; }
    std::shared_ptr<MultiTensor_<double>> get_contravariant(const int istate, const bool weight = false) const;

    const MatView shalf(const int iext, const int ist) const { return d_->shalf(to_denom_.at(iext), ist); }
    double phi(const int iext, const int ist, const size_t i) const { return d_->denom(to_denom_.at(iext), ist, i); }

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
