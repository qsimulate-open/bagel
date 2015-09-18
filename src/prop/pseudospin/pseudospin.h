//
// BAGEL - Parallel electron correlation program.
// Filename: pseudospin.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynolds2018@u.northwestern.edu>
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


#ifndef __SRC_PROP_PSEUDOSPIN_PSEUDOSPIN_H
#define __SRC_PROP_PSEUDOSPIN_PSEUDOSPIN_H

#include <src/wfn/geometry.h>
#include <src/util/math/zmatrix.h>
#include <src/util/kramers.h>
#include <src/wfn/rdm.h>
#include <src/ci/zfci/zharrison.h>

namespace bagel {

// Some operator that acts on a basis of pseudospin eigenstates
class Spin_Operator {
  protected:
    const int nspin_;
    std::shared_ptr<const ZMatrix> matrix_;
    std::complex<double> coeff_;

    // stevens_ == true:  order_ = k, index_ = q
    // stevens_ == false: order_ = 2, index_ from 0 to 8 denotes xx, yx, zx, xy...
    const bool stevens_;
    const int order_;
    const int index_;

  public:
    Spin_Operator(std::shared_ptr<const ZMatrix> _mat, const bool _stev, const int _ord, const int _ind);

    int nspin() const { return nspin_; }
    std::shared_ptr<const ZMatrix> matrix() const { return matrix_; }

    void set_coeff(const std::complex<double> in) { coeff_ = in; }
    std::complex<double> coeff() const { assert(coeff_ == coeff_); return coeff_; }  // Default value is NaN; triggers assertion if coeff_ has not been computed

    bool stevens() const { return stevens_; }
    int order() const { return order_; }
    int index() const { return index_; }

    std::string operator_name() const;
    std::string coeff_name() const;

    void print() const {
      const std::string label = operator_name() + " for S = " + (nspin_ == 1 ? "" : std::to_string(nspin_ / 2) + " ") + (nspin_ % 2 == 0 ? "" : "1/2");
      matrix_->print(label, 20);
    }
};


class Pseudospin {
  protected:
    int nspin_;
    int nspin1_;
    std::vector<double> ref_energy_;

    // These are the basic spin operators
    std::array<std::shared_ptr<ZMatrix>,3> spin_xyz_;
    std::shared_ptr<ZMatrix> spin_plus_;
    std::shared_ptr<ZMatrix> spin_minus_;

    // These are over eigenstates of the (ZFCI) Hamiltonian
    std::shared_ptr<ZMatrix> spinham_h_;
    std::array<std::shared_ptr<ZMatrix>,3> spinop_h_;

    // These are over pseudospin eigenstates
    std::shared_ptr<ZMatrix> spinham_s_;
    std::array<std::shared_ptr<ZMatrix>,3> spinop_s_;

    void update_spin_matrices(VectorB spinvals);

  public:
    Pseudospin(const int nspin);

    std::vector<Spin_Operator> build_extended_stevens_operators(const std::vector<int> ranks) const;
    std::vector<Spin_Operator> build_2ndorder_zfs_operators() const;

    void compute_numerical_hamiltonian(const ZHarrison& zfci, std::shared_ptr<const RelCoeff_Block> active_coeff);
    std::vector<Spin_Operator> extract_hamiltonian_parameters(const bool real, const std::vector<Spin_Operator> param);

    std::shared_ptr<ZMatrix> spin_xyz(const int i) const { return spin_xyz_[i]; }
    std::shared_ptr<ZMatrix> spin_plus() const { return spin_plus_; }
    std::shared_ptr<ZMatrix> spin_minus() const { return spin_minus_; }

    static std::shared_ptr<ZMatrix> compute_Dtensor(const std::vector<Spin_Operator> input);
};

}

#endif
