//
// BAGEL - Parallel electron correlation program.
// Filename: shell.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_MOLECULE_SHELL_H
#define __SRC_MOLECULE_SHELL_H

#include <array>
#include <src/math/matrix.h>

namespace bagel {

class Shell {

  protected:
    const bool spherical_;

    std::array<double,3> position_;
    const int angular_number_;
    const std::vector<double> exponents_;     // length of primitive basis function
    const std::vector<std::vector<double>> contractions_;  // length of contracted basis function
    const std::vector<std::pair<int, int>> contraction_ranges_;

    const bool dummy_;
    std::vector<int> contraction_upper_;
    std::vector<int> contraction_lower_;

    int nbasis_;

    // whether a relativistic part is initialized
    bool relativistic_;

    // protected members for relativistic calculations
    std::array<std::shared_ptr<const Matrix>,3> small_;
    std::shared_ptr<const Shell> aux_increment_;
    std::shared_ptr<const Shell> aux_decrement_;
    std::shared_ptr<const Matrix> overlap_compute_() const;
    std::array<std::shared_ptr<const Matrix>,3> moment_compute_(const std::shared_ptr<const Matrix> overlap) const;

  public:
    Shell(const bool spherical, const std::array<double,3>& position, int angular_num, const std::vector<double>& exponents,
          const std::vector<std::vector<double>>& contraction, const std::vector<std::pair<int, int>>& cont_range);
    // default constructor for adding null basis
    Shell(const bool sph);

    bool dummy() const { return dummy_; };
    bool spherical() const { return spherical_; };
    int num_primitive() const { return exponents_.size(); };
    int num_contracted() const { return contractions_.size(); };

    double position(const int i) const { return position_[i]; };
    const std::array<double,3>& position() const { return position_; };
    int angular_number() const { return angular_number_; };
    double exponents(const int i) const { return exponents_[i]; };
    const std::vector<double>& exponents() const { return exponents_; };
    const double* exponents_pointer() const { return &(exponents_[0]); };
    const std::vector<double> contractions(const int i) const { return contractions_[i]; };
    const std::vector<std::vector<double>>& contractions() const { return contractions_; };
    const std::pair<int, int>& contraction_ranges(const int i) const { return contraction_ranges_[i]; };
    const std::vector<std::pair<int, int>>& contraction_ranges() const { return contraction_ranges_; };

    const std::vector<int>& contraction_upper() const { return contraction_upper_; };
    const std::vector<int>& contraction_lower() const { return contraction_lower_; };

    const std::string show() const;
    int nbasis() const { return nbasis_; };

    std::shared_ptr<const Shell> move_atom(const std::array<double,3>&) const;
    std::shared_ptr<const Shell> move_atom(const double*) const;

    bool operator==(const Shell& o) const;

    // generates a shell that satisfy kinetic balance at the primitive level.
    template<int inc>
    std::shared_ptr<const Shell> kinetic_balance_uncont() const;

    std::shared_ptr<const Shell> cartesian_shell() const;

    std::vector<std::shared_ptr<const Shell>> split_if_possible(const size_t batchsize) const;

    void init_relativistic();

    // Relativistic
    bool relativistic() const { return relativistic_; }
    const std::shared_ptr<const Matrix> small(const int i) const { assert(relativistic_); return small_[i]; }
    const std::shared_ptr<const Shell> aux_increment() const { assert(relativistic_); return aux_increment_; }
    const std::shared_ptr<const Shell> aux_decrement() const { assert(relativistic_); return aux_decrement_; }
    int nbasis_aux_increment() const { return aux_increment_ ? aux_increment_->nbasis() : 0; }
    int nbasis_aux_decrement() const { return aux_decrement_ ? aux_decrement_->nbasis() : 0; }

    // DFT grid
    void compute_grid_value(double*, double*, double*, double*, const double& x, const double& y, const double& z) const;
    void compute_grid_value_deriv2(double*, double*, double*, double*, double*, double*, const double& x, const double& y, const double& z) const;
};

}

#endif

