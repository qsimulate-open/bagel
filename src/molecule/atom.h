//
// BAGEL - Parallel electron correlation program.
// Filename: atom.h
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


#ifndef __SRC_MOLECULE_ATOM_H
#define __SRC_MOLECULE_ATOM_H

#include <src/molecule/shell.h>
#include <src/input/input.h>

namespace bagel {

class Atom {
  protected:

    bool spherical_;

    std::string name_;
    std::array<double,3> position_;
    std::vector<std::shared_ptr<const Shell>> shells_;
    int atom_number_;
    double atom_charge_;
    double atom_exponent_;
    int nbasis_;
    int lmax_;

    // basis set
    std::string basis_;

    // effective core potential (0 = zeta, 1 = coef)
    std::array<double,2> ecp_;

    // This function sets shell_ and lmax_
    // in : a vector of an angular label, exponents, and coefficients.
    void construct_shells(std::vector<std::tuple<std::string, std::vector<double>, std::vector<std::vector<double>>>> in);
    // if needed and possible, we split shells whose nbasis are bigger than batchsize
    void split_shells(const size_t batchsize);

    void basis_init(std::shared_ptr<const PTree>);
    void common_init();

  public:
    Atom(std::shared_ptr<const PTree> inp, const bool spherical, const bool angstrom, const std::pair<std::string, std::shared_ptr<const PTree>> defbas, std::shared_ptr<const PTree> elem, const bool aux= false);

    Atom(const bool spherical, const std::string name, const std::array<double,3>& position, const std::string bas, const std::pair<std::string, std::shared_ptr<const PTree>> json, std::shared_ptr<const PTree> elem);
    Atom(const bool spherical, const std::string name, const std::array<double,3>& position, const double charge);
    Atom(const bool spherical, const std::string name, const std::array<double,3>& position,
         const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>);
    Atom(const std::string name, const std::string bas, const std::vector<std::shared_ptr<const Shell>> shell);

    Atom(const Atom&, const bool spherical, const std::string bas, const std::pair<std::string, std::shared_ptr<const PTree>> defbas, std::shared_ptr<const PTree> elem);
    Atom(const Atom&, const std::array<double,3>& displ);

    const std::string name() const { return name_; }
    int atom_number() const { return atom_number_;}
    double atom_charge() const { return atom_charge_;}

    double atom_exponent() const { return atom_exponent_; }
    bool finite_nucleus() const { return atom_exponent_ != 0.0; }

    double radius() const;
    double cov_radius() const;
    const std::array<double,3>& position() const { return position_; }
    double position(const unsigned int i) const { return position_[i]; }
    const std::vector<std::shared_ptr<const Shell>>& shells() const { return shells_; }
    int nshell() const { return shells_.size(); }

    bool dummy() const { return atom_number_ == 0; }

    int nbasis() const { return nbasis_; }
    int lmax() const { return lmax_; }
    bool spherical() const { return spherical_; }

    std::string basis() const { return basis_; }

    void print_basis() const;
    void print() const;

    bool operator==(const Atom&) const;

    // effective core potential
    bool ecp() const { return ecp(0) != 0.0 || ecp(1) != 0.0; }
    double ecp(const int i) const { return ecp_[i]; }

    // distance between this and other
    double distance(const std::array<double,3>& o) const;
    double distance(const std::shared_ptr<const Atom> o) const { return distance(o->position()); }
    // displacement vector from this to other
    std::array<double,3> displ(const std::shared_ptr<const Atom> o) const;
    // angle between two atoms (A,B)  and this (O) deg(A-O-B)
    double angle(const std::shared_ptr<const Atom>, const std::shared_ptr<const Atom>) const;
    // dihedral angle for A-this-O-B
    double dihedral_angle(const std::shared_ptr<const Atom>, const std::shared_ptr<const Atom>, const std::shared_ptr<const Atom>) const;

    // inititalize relativistic calculation
    std::shared_ptr<const Atom> relativistic() const;
};

}

#endif

