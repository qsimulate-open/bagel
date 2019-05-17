//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: molecule.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_MOLECULE_MOLECULE_H
#define __SRC_MOLECULE_MOLECULE_H

#include <iostream>
#include <src/molecule/atom.h>
#include <src/util/math/xyzfile.h>
#include <src/util/serialization.h>
#include <src/opt/constraint.h>

namespace bagel {

class Molecule {
  protected:
    bool spherical_;

    bool aux_merged_;

    // Some shared info for basis sets.
    int nbasis_;
    int nele_;
    int nfrc_;
    int naux_;
    int lmax_;
    int aux_lmax_;

    // these two offsets are in principle redundant information (can be derived from Shells);
    std::vector<std::vector<int>> offsets_;
    std::vector<std::vector<int>> aux_offsets_;

    std::string basisfile_;
    std::string auxfile_;

    // Atoms, which contains basis-set info also.
    std::vector<std::shared_ptr<const Atom>> atoms_;
    std::vector<std::shared_ptr<const Atom>> aux_atoms_;

    // Nuclear repulsion energy.
    double nuclear_repulsion_;

    // external electric field
    std::array<double,3> external_;

    // external magnetic field
    std::array<double,3> magnetic_field_;

    // Computes the nuclear repulsion energy.
    bool skip_self_interaction_;
    double compute_nuclear_repulsion();

    // Constructor helpers
    void common_init1();

    // Absorbing potential
    std::array<double,3> cap_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & spherical_ & aux_merged_ & nbasis_ & nele_ & nfrc_ & naux_ & lmax_ & aux_lmax_ & offsets_ & aux_offsets_ & basisfile_ & auxfile_
         & atoms_ & aux_atoms_ & nuclear_repulsion_ & external_ & magnetic_field_ & skip_self_interaction_;
    }

  public:
    Molecule() : external_{{0.0,0.0,0.0}}, magnetic_field_{{0.0,0.0,0.0}} {}
    Molecule(const std::vector<std::shared_ptr<const Atom>> a, const std::vector<std::shared_ptr<const Atom>> b)
      : atoms_(a), aux_atoms_(b), external_{{0.0,0.0,0.0}}, magnetic_field_{{0.0,0.0,0.0}} { }
    Molecule(const Molecule& o, std::shared_ptr<const Matrix> disp, const bool rotate = true);
    virtual ~Molecule() { }

    // Returns shared pointers of Atom objects, which contains basis-set info.
    const std::vector<std::shared_ptr<const Atom>>& atoms() const { return atoms_; }
    const std::vector<std::shared_ptr<const Atom>>& aux_atoms() const { return aux_atoms_; }
    std::shared_ptr<const Atom> atoms(const unsigned int i) const { return atoms_[i]; }

    // Returns constants and private members
    bool spherical() const { return spherical_; }
    size_t nbasis() const { return nbasis_; }
    size_t nele() const { return nele_; }
    size_t nfrc() const { return nfrc_; }
    size_t naux() const { return naux_; }
    int lmax() const { return lmax_; }
    int aux_lmax() const { return aux_lmax_; }
    int natom() const { return atoms_.size(); }
    const std::string basisfile() const { return basisfile_; }
    const std::string auxfile() const { return auxfile_; }
    virtual double nuclear_repulsion() const { return nuclear_repulsion_; }
    bool skip_self_interaction() { return skip_self_interaction_; }

    // The position of the specific function in the basis set.
    const std::vector<std::vector<int>>& offsets() const { return offsets_; }
    const std::vector<std::vector<int>>& aux_offsets() const { return aux_offsets_; }
    const std::vector<int>& offset(const unsigned int i) const { return offsets_.at(i); }
    const std::vector<int>& aux_offset(const unsigned int i) const { return aux_offsets_.at(i); }

    int num_count_ncore_only() const; // also set nfrc_

    void print_atoms() const;

    bool operator==(const Molecule& o) const;

    std::array<double,3> charge_center() const;
    std::array<double,6> quadrupole() const;

    // finite nucleus
    bool has_finite_nucleus() const;

    // external electric field
    bool external() const { return external(0) != 0.0 || external(1) != 0.0 || external(2) != 0.0; }
    double external(const int i) const { return external_[i]; }

    // external magnetic field
    bool nonzero_magnetic_field() const { return magnetic_field(0) != 0.0 || magnetic_field(1) != 0.0 || magnetic_field(2) != 0.0; }
    std::array<double,3> magnetic_field() const { return magnetic_field_; }
    double magnetic_field(const int i) const { return magnetic_field_[i]; }

    std::shared_ptr<const XYZFile> xyz() const;

    std::array<double,3> cap() const { return cap_; }

    // In R12 methods, we need to construct a union of OBS and CABS.
    // Currently, this is done by creating another object and merge OBS and CABS into atoms_.
    // After this, compute_nuclear_repulsion() should not be called.
    // Not undo-able.
    void merge_obs_aux();

    // transformation matrices for the internal coordinate for geometry optimization
    // ninternal runs fast (and cartsize slower) (weighted Wilson B)
    std::array<std::shared_ptr<const Matrix>,3> compute_internal_coordinate(
        std::shared_ptr<const Matrix> prev = nullptr,
        std::vector<std::shared_ptr<const OptExpBonds>> explicit_bond = std::vector<std::shared_ptr<const OptExpBonds>>(),
        bool negative_hessian = false,
        bool verbose = true) const;
    // driver for compute B matrix for redundant coordinate (original Wilson B)
    std::tuple<std::vector<std::vector<int>>,std::array<std::shared_ptr<const Matrix>,5>> compute_redundant_coordinate(std::vector<std::vector<int>> prev) const;

    // Split up the atoms into several Molecule objects
    // To limit the memory requirement of integral evaluation
    const std::vector<std::shared_ptr<const Molecule>> split_atoms(const int max_atoms) const;

    std::shared_ptr<Molecule> uncontract() const;

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Molecule)

namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<Molecule, T>::value>::type> {
    typedef Molecule type;
  };
}

#endif
