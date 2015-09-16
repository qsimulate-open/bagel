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
    int nspin_;
    std::shared_ptr<const ZMatrix> matrix_;
    std::string coeff_name_;

  public:
    Spin_Operator(std::shared_ptr<const ZMatrix> mat, const std::string name)
     : matrix_(mat), coeff_name_(name) {
      assert(matrix_->ndim() == matrix_->mdim());
      nspin_ = matrix_->ndim();
    }

    int nspin() const { return nspin_; }
    std::shared_ptr<const ZMatrix> matrix() const { return matrix_; }
    std::string coeff_name() const { return coeff_name_; }
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
    void extract_hamiltonian_parameters(const bool real, const std::vector<Spin_Operator> param);

    std::shared_ptr<ZMatrix> spin_xyz(const int i) const { return spin_xyz_[i]; }
    std::shared_ptr<ZMatrix> spin_plus() const { return spin_plus_; }
    std::shared_ptr<ZMatrix> spin_minus() const { return spin_minus_; }
};

}

#endif
