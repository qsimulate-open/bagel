//
// BAGEL - Parallel electron correlation program.
// Filename: pgeometry.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __src_pscf_pgeometry_h
#define __src_pscf_pgeometry_h

#include <string>
#include <src/scf/geometry.h>

namespace bagel {

class PGeometry : public Geometry {
  protected:
    int L_; // Namur cutoff L
    int S_; // Namur cutoff S
    int K_; // number of K points (= number of unit cell in the *half* first Brillouin zone)
    double A_; // fundamental vector

    // Computes the nuclear repulsion energy per unit cell.
    double pnuclear_repulsion_;
    double compute_pnuclear_repulsion() const;

  public:
    PGeometry(const std::string, const int);
    ~PGeometry() {};

    // Some constants for periodic calculations.
    int L() const { return L_; };
    int S() const { return S_; };
    int K() const { return K_; };
    double A() const { return A_; };

    // Returns nuclear repulsion energies.
    double nuclear_repulsion() const { return pnuclear_repulsion_; };
};

}

#endif

