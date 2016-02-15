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
class Stevens_Operator {
  protected:
    const int nspin_;
    std::shared_ptr<const ZMatrix> matrix_;
    double coeff_;

    // order_ = k, index_ = q
    const int order_;
    const int index_;

  public:
    Stevens_Operator(std::shared_ptr<const ZMatrix> _mat, const int _ord, const int _ind);

    int nspin() const { return nspin_; }
    std::shared_ptr<const ZMatrix> matrix() const { return matrix_; }

    void set_coeff(const double in) { coeff_ = in; }
    double coeff() const { assert(coeff_ == coeff_); return coeff_; }  // Default value is NaN; triggers assertion if coeff_ has not been computed

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
    std::vector<Stevens_Operator> ESO_;
    std::shared_ptr<const Matrix> dtensor_;
    double Dval_;
    double Eval_;

    std::shared_ptr<const PTree> idata_;

    // These are the basic spin operators
    std::array<std::shared_ptr<ZMatrix>,3> spin_xyz_;
    std::shared_ptr<ZMatrix> spin_plus_;
    std::shared_ptr<ZMatrix> spin_minus_;

    // Columns of this matrix give the axes along which we define Sx, Sy, and Sz
    // Normally, they are eigenvectors of the Abragam-Bleaney tensor (G)
    std::shared_ptr<const Matrix> spin_axes_;

    // These matrices run over eigenstates of the (ZFCI) Hamiltonian; those over pseudospin eigenstates are not saved as members
    std::shared_ptr<ZMatrix> spinham_h_;
    std::shared_ptr<ZMatrix> trev_h_;
    std::array<std::shared_ptr<ZMatrix>,3> zfci_mu_;
    std::array<std::shared_ptr<ZMatrix>,3> zfci_spin_;
    std::array<std::shared_ptr<ZMatrix>,3> zfci_orbang_;

    // Same as the above three, but rotated from geometric axes to spin_axes_
    std::array<std::shared_ptr<ZMatrix>, 3> zfci2_mu_;
    std::array<std::shared_ptr<ZMatrix>, 3> zfci2_spin_;
    std::array<std::shared_ptr<ZMatrix>, 3> zfci2_orbang_;

    void update_spin_matrices(VectorB spinvals);
    std::shared_ptr<const Matrix> read_axes(std::shared_ptr<const Matrix> default_axes) const;

  public:
    Pseudospin(const int nspin, std::shared_ptr<const PTree> idata);

    void compute(const ZHarrison& zfci);

    // return symbolic spin matrices
    std::shared_ptr<ZMatrix> spin_xyz(const int i) const { return spin_xyz_[i]; }
    std::shared_ptr<ZMatrix> spin_plus() const { return spin_plus_; }
    std::shared_ptr<ZMatrix> spin_minus() const { return spin_minus_; }

    // return calculated ZFS parameters
    std::vector<Stevens_Operator> ESO() const { return ESO_; }
    Stevens_Operator ESO(const int i) const { return ESO_[i]; }
    std::shared_ptr<const Matrix> dtensor() const { return dtensor_; }
    double Dval() const { return Dval_; }
    double Eval() const { return Eval_; }

    // setup functions
    std::vector<Stevens_Operator> build_extended_stevens_operators(const std::vector<int> ranks) const;
    void compute_numerical_hamiltonian(const ZHarrison& zfci, std::shared_ptr<const RelCoeff_Block> active_coeff);
    std::shared_ptr<const Matrix> identify_magnetic_axes() const;

    // to extract D-tensor
    std::shared_ptr<const ZMatrix> compute_spin_eigenvalues() const;
    std::vector<Stevens_Operator> extract_hamiltonian_parameters(const std::vector<Stevens_Operator> param, std::shared_ptr<const ZMatrix> spinham_s) const;
    std::tuple<std::shared_ptr<const Matrix>, double, double> compute_Dtensor(const std::vector<Stevens_Operator> input) const;

    // To verify that a matrix in pseudospin basis has specified time-reversal symmetry
    // The bool parameters tell us if the matrix should be symmetric or antisymmetric under Hermitian conjugation and time-reversal, respectively
    bool is_t_symmetric(const ZMatrix& in, const bool hermitian = true, const bool t_symmetric = true, const double thresh = 1.0e-6) const;

};

}

#endif
