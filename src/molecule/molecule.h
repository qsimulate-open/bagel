//
// BAGEL - Parallel electron correlation program.
// Filename: molecule.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_MOLECULE_MOLECULE_H
#define __SRC_MOLECULE_MOLECULE_H

#include <iostream>
#include <src/molecule/atom.h>
#include <src/molecule/petite.h>
#include <src/math/xyzfile.h>

namespace bagel {

class Molecule {
  protected:
    // Atoms, which contains basis-set info also.
    std::vector<std::shared_ptr<const Atom>> atoms_;
    std::vector<std::shared_ptr<const Atom>> aux_atoms_;

    // Nuclear repulsion energy.
    double nuclear_repulsion_;

    // Symmetry can be used for molecular calculation.
    std::string symmetry_;
    std::shared_ptr<Petite> plist_;
    int nirrep_;

    // external field
    std::array<double,3> external_;

    // Computes the nuclear repulsion energy.
    double compute_nuclear_repulsion();

  public:
    Molecule() : symmetry_("c1"), nirrep_(1) {}
    Molecule(const std::vector<std::shared_ptr<const Atom>> a, const std::vector<std::shared_ptr<const Atom>> b) : atoms_(a), aux_atoms_(b), symmetry_("c1"), nirrep_(1) { }

    // Returns shared pointers of Atom objects, which contains basis-set info.
    const std::vector<std::shared_ptr<const Atom>>& atoms() const { return atoms_; }
    const std::vector<std::shared_ptr<const Atom>>& aux_atoms() const { return aux_atoms_; }
    std::shared_ptr<const Atom> atoms(const unsigned int i) const { return atoms_[i]; }

    // Returns constants and private members
    int natom() const { return atoms_.size(); }
    virtual double nuclear_repulsion() const { return nuclear_repulsion_; }
    const std::string symmetry() const { return symmetry_; }

    void print_atoms() const;

    std::array<double,3> charge_center() const;
    std::array<double,6> quadrupole() const;

    // finite nucleus
    bool has_finite_nucleus() const;

    // external field
    bool external() const { return external(0) != 0.0 || external(1) != 0.0 || external(2) != 0.0; }
    double external(const int i) const { return external_[i]; }

    virtual size_t nbasis() const { return std::accumulate(atoms_.begin(), atoms_.end(), 0,
                                    [](const int& i, const std::shared_ptr<const Atom>& j) { return i+j->nbasis(); }); }
    virtual size_t naux() const { return std::accumulate(aux_atoms_.begin(), aux_atoms_.end(), 0,
                                    [](const int& i, const std::shared_ptr<const Atom>& j) { return i+j->nbasis(); }); }

    std::shared_ptr<const XYZFile> xyz() const;

};

}

#endif
