//
// BAGEL - Parallel electron correlation program.
// Filename: atomgroup.h
// Copyright (C) 2015 Toru Shiozaki
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


#ifndef __SRC_PERIODIC_ATOMGROUP_H
#define __SRC_PERIODIC_ATOMGROUP_H

#include <src/molecule/atom.h>

namespace bagel {

class AtomGroup {
  protected:
    std::vector<std::shared_ptr<const Atom>> atoms_;
    std::vector<int> order_in_geom_, ishell_;
    std::array<double, 3> position_;
    int nbasis_, natom_;

    void init() {
      position_ = {{0.0, 0.0, 0.0}};
      nbasis_ = 0;
      double sum = 0.0;
      for (auto& atom : atoms_) {
        position_[0] += atom->atom_charge() * atom->position(0);
        position_[1] += atom->atom_charge() * atom->position(1);
        position_[2] += atom->atom_charge() * atom->position(2);
        sum += atom->atom_charge();
        nbasis_ += atom->nbasis();
      }
      position_[0] /= sum;
      position_[1] /= sum;
      position_[2] /= sum;
    }

  public:
    AtomGroup(std::vector<std::shared_ptr<const Atom>> satom, std::vector<int> order, std::vector<int> ishell)
     : atoms_(satom), order_in_geom_(order), ishell_(ishell) {

      init();
    }
    ~AtomGroup() { }

    std::array<double, 3> position() const { return position_; }
    double position(const int i) const { return position_[i]; }

    std::vector<std::shared_ptr<const Atom>> atoms() const { return atoms_; }
    std::shared_ptr<const Atom> atoms(const int i) const { return atoms_[i]; }
    std::vector<int> order_in_geom()  const { return order_in_geom_; }
    int order_in_geom(const int i)  const { return order_in_geom_[i]; }
    std::vector<int> ishell()  const { return ishell_; }
    int ishell(const int i)  const { return ishell_[i]; }
    int nbasis() const { return nbasis_; }
    int natom() const { return atoms_.size(); }
};

}
#endif
