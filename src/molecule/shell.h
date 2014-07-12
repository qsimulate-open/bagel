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

#include <src/math/matrix.h>
#include <src/math/zmatrix.h>
#include <src/molecule/shell_base.h>
#include <src/util/serialization.h>

namespace bagel {

class Shell : public Shell_base {

  protected:
    std::vector<double> exponents_;     // length of primitive basis function
    std::vector<std::vector<double>> contractions_;  // length of contracted basis function
    std::vector<std::pair<int, int>> contraction_ranges_;

    bool dummy_;
    std::vector<int> contraction_upper_;
    std::vector<int> contraction_lower_;

    int nbasis_;

    // whether a relativistic part is initialized
    bool relativistic_;

    // protected members for relativistic calculations
    std::array<std::shared_ptr<const Matrix>,3>  small_;
    std::shared_ptr<const Shell> aux_increment_;
    std::shared_ptr<const Shell> aux_decrement_;

    // TODO Refactor - These next few are essentially the same as above, but for London integrals only
    std::shared_ptr<const Shell> aux_same_;
    std::array<std::shared_ptr<const ZMatrix>,3> zsmall_;
    std::array<std::shared_ptr<const ZMatrix>,3> zsmallc_;

    // whether a london phase factor is being used
    bool magnetism_;
    std::array<double,3> vector_potential_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << spherical_ << position_ << vector_potential_ << angular_number_ << exponents_ << contractions_ << contraction_ranges_
         << dummy_ << contraction_upper_ << contraction_lower_ << nbasis_ << relativistic_;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> spherical_ >> position_ >> vector_potential_ >> angular_number_ >> exponents_ >> contractions_ >> contraction_ranges_
         >> dummy_ >> contraction_upper_ >> contraction_lower_ >> nbasis_ >> relativistic_;
      if (relativistic_)
        init_relativistic();
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }

  public:
    Shell() { }
    Shell(const bool spherical, const std::array<double,3>& position, int angular_num, const std::vector<double>& exponents,
          const std::vector<std::vector<double>>& contraction, const std::vector<std::pair<int, int>>& cont_range);
    // default constructor for adding null basis
    Shell(const bool sph);

    bool dummy() const { return dummy_; };
    int num_primitive() const { return exponents_.size(); };
    int num_contracted() const { return contractions_.size(); };

    double exponents(const int i) const { return exponents_[i]; };
    const std::vector<double>& exponents() const { return exponents_; };
    const double* exponents_pointer() const { return &(exponents_[0]); };
    const std::vector<double> contractions(const int i) const { return contractions_[i]; };
    const std::vector<std::vector<double>>& contractions() const { return contractions_; };
    const std::pair<int, int>& contraction_ranges(const int i) const { return contraction_ranges_[i]; };
    const std::vector<std::pair<int, int>>& contraction_ranges() const { return contraction_ranges_; };

    const std::vector<int>& contraction_upper() const { return contraction_upper_; };
    const std::vector<int>& contraction_lower() const { return contraction_lower_; };

    int nbasis() const { return nbasis_; };
    std::string show() const override;

    std::shared_ptr<const Shell> move_atom(const std::array<double,3>&) const;
    std::shared_ptr<const Shell> move_atom(const double*) const;

    bool operator==(const Shell& o) const;

    // generates a shell that satisfy kinetic balance at the primitive level.
    template<int inc>
    std::shared_ptr<const Shell> kinetic_balance_uncont() const;

    std::shared_ptr<const Shell> cartesian_shell() const;

    std::vector<std::shared_ptr<const Shell>> split_if_possible(const size_t batchsize) const;

    // magnetism
    bool magnetism() const { return magnetism_; }
    double vector_potential(const unsigned int i) const { return vector_potential_[i]; }
    const std::array<double,3>& vector_potential() const { return vector_potential_; };
    void add_phase(const std::array<double,3>& phase_input);

    void init_relativistic();
    void init_relativistic_london(const std::array<double,3> magnetic_field, bool london);

    // Relativistic
    bool relativistic() const { return relativistic_; }
    const std::shared_ptr<const Matrix>  small(const int i)  const { assert(relativistic_); return  small_[i]; }
    const std::shared_ptr<const ZMatrix> zsmall(const int i) const { assert(relativistic_); return zsmall_[i]; }
    const std::shared_ptr<const ZMatrix> zsmallc(const int i) const { assert(relativistic_); return zsmallc_[i]; }
    const std::shared_ptr<const Shell> aux_increment() const { assert(relativistic_); return aux_increment_; }
    const std::shared_ptr<const Shell> aux_decrement() const { assert(relativistic_); return aux_decrement_; }
    const std::shared_ptr<const Shell> aux_same() const { assert(relativistic_); return aux_same_; }
    int nbasis_aux_increment() const { return aux_increment_ ? aux_increment_->nbasis() : 0; }
    int nbasis_aux_decrement() const { return aux_decrement_ ? aux_decrement_->nbasis() : 0; }
    int nbasis_aux_same() const { return aux_same_ ? aux_same_->nbasis() : 0; }

    // DFT grid
    void compute_grid_value(double*, double*, double*, double*, const double& x, const double& y, const double& z) const;
    void compute_grid_value_deriv2(double*, double*, double*, double*, double*, double*, const double& x, const double& y, const double& z) const;
};

}

#endif

