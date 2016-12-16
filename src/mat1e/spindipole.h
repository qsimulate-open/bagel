//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: spindipole.h
// Copyright (C) 2015 Toru Shiozaki
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


#ifndef __SRC_MAT1E_SPINDIPOLE_H
#define __SRC_MAT1E_SPINDIPOLE_H

#include <src/mat1e/matrix1earray.h>

namespace bagel {

class SpinDipole : public Matrix1eArray<6> {
  protected:
    std::shared_ptr<const Atom> atom_;
    int s_;

    void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int, std::shared_ptr<const Molecule>) override;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Matrix1e>(*this);
    }

  public:
    SpinDipole() { }
    SpinDipole(std::shared_ptr<const Molecule>, std::shared_ptr<const Atom>, const int s);

    void print(const std::string name, const int len = 10) const override { Matrix1eArray<6>::print(name.empty() ? "Spin-Dipole" : name, len); }

};

}

#endif

