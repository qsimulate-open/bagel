//
// BAGEL - Parallel electron correlation program.
// Filename: lattice.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __BAGEL_SRC_PERIODIC_LATTICE_H
#define __BAGEL_SRC_PERIODIC_LATTICE_H

#include <src/wfn/geometry.h>

namespace bagel {

class Lattice {
  protected:
    int ndim_;
    int ncell_; // tmp
    int num_lattice_pts_;
    std::shared_ptr<const Geometry> primitive_cell_;
    std::vector<std::array<double, 3>> lattice_vectors_;

    double nuclear_repulsion_;
    double compute_nuclear_repulsion() const;

  public:
    Lattice() { }
    Lattice(const std::shared_ptr<const Geometry> g);
    virtual ~Lattice() { }

    int ndim() const { return ndim_; }
    int ncell() const {return ncell_; }
    int num_lattice_pts() const { return num_lattice_pts_; }
    std::shared_ptr<const Geometry> primitive_cell() const { return primitive_cell_; }
    std::vector<std::array<double, 3>> lattice_vectors() const { return lattice_vectors_; }
    std::array<double, 3> lattice_vectors(const int i) const { return lattice_vectors_[i]; }

    void init();
    double nuclear_repulsion() const { return nuclear_repulsion_; };
    void print_primitive_vectors() const;
    void print_lattice_coordinates() const; // write .XYZ file
};

}

#endif

