//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: atom.h
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


#ifndef __SRC_MOLECULE_ATOM_H
#define __SRC_MOLECULE_ATOM_H

#include <src/molecule/shell.h>
#include <src/molecule/ecp.h>
#include <src/molecule/soecp.h>
#include <src/util/input/input.h>

namespace bagel {

class Atom {
  protected:
    bool spherical_;

    std::string name_;
    std::array<double,3> position_;
    std::array<double,3> vector_potential_;
    std::vector<std::shared_ptr<const Shell>> shells_;
    bool use_ecp_basis_;
    std::shared_ptr<const ECP> ecp_parameters_;
    std::shared_ptr<const SOECP> so_parameters_;
    int atom_number_;
    double atom_charge_;
    double atom_exponent_;
    double mass_;
    int nbasis_;
    int lmax_;

    // basis set
    std::string basis_;

    // This function sets shell_ and lmax_
    // in : a vector of an angular label, exponents, and coefficients.
    void construct_shells(std::vector<std::tuple<std::string, std::vector<double>, std::vector<std::vector<double>>>> in);
    // in : angular momentum (l), exponents (zeta_kl), coefficients (A_kl), powers of r (n_kl)
    void construct_shells_ECP(const int ncore, std::vector<std::tuple<std::string, std::vector<double>,
                              std::vector<double>, std::vector<int>>> in);
    void construct_shells_SOECP(std::vector<std::tuple<std::string, std::vector<double>,
                                std::vector<double>, std::vector<int>>> in);

    // if needed and possible, we split shells whose nbasis are bigger than batchsize
    void split_shells(const size_t batchsize);

    void basis_init(std::shared_ptr<const PTree>);
    void basis_init_ECP(std::shared_ptr<const PTree>);
    void common_init();

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & spherical_ & name_ & position_ & vector_potential_ & shells_ & use_ecp_basis_ & ecp_parameters_ & so_parameters_
         & atom_number_ & atom_charge_ & atom_exponent_ & mass_ & nbasis_ & lmax_ & basis_;
    }

  public:
    Atom() { }
    Atom(std::shared_ptr<const PTree> inp, const bool spherical, const bool angstrom, const std::pair<std::string, std::shared_ptr<const PTree>> defbas,
         std::shared_ptr<const PTree> elem, const bool aux=false, const bool ecp=false, const bool default_finite=false);

    Atom(const bool spherical, const std::string name, const std::array<double,3>& position, const std::string bas,
         const std::pair<std::string, std::shared_ptr<const PTree>> json, std::shared_ptr<const PTree> elem);
    Atom(const bool spherical, const std::string name, const std::array<double,3>& position, const double charge);
    Atom(const bool spherical, const std::string name, const std::array<double,3>& position,
         const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>, const std::string bas = "custom_basis");
    Atom(const std::string name, const std::string bas, const std::vector<std::shared_ptr<const Shell>> shell);
    Atom(const std::string name, const std::string bas, const std::vector<std::shared_ptr<const Shell>> shell,
                                                        const std::vector<std::shared_ptr<const Shell_ECP>> shell_ECP, const int ncore, const int maxl);
    Atom(const std::string name, const std::string bas, const std::vector<std::shared_ptr<const Shell>> shell,
                                                        const std::shared_ptr<const ECP> ecp_param);

    Atom(const Atom&, const bool spherical, const std::string bas, const std::pair<std::string, std::shared_ptr<const PTree>> defbas,
         std::shared_ptr<const PTree> elem);
    Atom(const Atom&, const std::array<double,3>& displ);

    const std::string name() const { return name_; }
    int atom_number() const { return atom_number_;}
    double atom_charge() const { return atom_charge_;}
    double atom_exponent() const { return atom_exponent_; }
    double mass() const { return mass_; }
    bool finite_nucleus() const { return atom_exponent_ != 0.0; }

    double radius() const;
    double cov_radius() const;
    const std::array<double,3>& position() const { return position_; }
    double position(const unsigned int i) const { return position_[i]; }
    const std::vector<std::shared_ptr<const Shell>>& shells() const { return shells_; }
    int nshell() const { return shells_.size(); }

    bool use_ecp_basis() const { return use_ecp_basis_; }
    const std::shared_ptr<const ECP>& ecp_parameters() const { return ecp_parameters_; }
    const std::shared_ptr<const SOECP>& so_parameters() const { return so_parameters_; }

    bool dummy() const { return atom_number_ == 0; }

    int nbasis() const { return nbasis_; }
    int lmax() const { return lmax_; }
    bool spherical() const { return spherical_; }

    std::string basis() const { return basis_; }

    void print_basis() const;
    void print() const;

    bool operator==(const Atom&) const;

    // distance between this and other
    double distance(const std::array<double,3>& o) const;
    double distance(const std::shared_ptr<const Atom> o) const { return distance(o->position()); }
    // displacement vector from this to other
    std::array<double,3> displ(const std::shared_ptr<const Atom> o) const;
    // angle between two atoms (A,B)  and this (O) deg(A-O-B)
    double angle(const std::shared_ptr<const Atom>, const std::shared_ptr<const Atom>) const;
    // dihedral angle for A-this-O-B
    double dihedral_angle(const std::shared_ptr<const Atom>, const std::shared_ptr<const Atom>, const std::shared_ptr<const Atom>) const;

    double vector_potential(const unsigned int i) const { assert(false);  return vector_potential_[i]; }
    const std::array<double,3>& vector_potential() const { assert(false);  return vector_potential_; };

    // inititalize relativistic calculation
    std::shared_ptr<const Atom> relativistic() const;
    std::shared_ptr<const Atom> relativistic(const std::array<double,3>& magnetic_field, bool london) const;

    // initialize magnetic field calculations
    std::shared_ptr<const Atom> apply_magnetic_field(const std::array<double,3>& field, const bool london) const;

    std::shared_ptr<const Atom> uncontract() const;

    void reset_shells(std::vector<std::shared_ptr<const Shell>>);

};

}

#endif

