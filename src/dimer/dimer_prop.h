
// BAGEL - Parallel electron correlation program.
// Filename: dimer_prop.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __BAGEL_DIMER_PROP_H
#define __BAGEL_DIMER_PROP_H

#include <array>
#include <src/fci/properties.h>

namespace bagel {

class DimerDipole : public CIDipole {
  protected:
    std::array<std::shared_ptr<const Matrix>, 3> dipole_moA_;
    std::array<std::shared_ptr<const Matrix>, 3> dipole_moB_;

    std::array<std::shared_ptr<const Matrix>, 3> cross_dipole_;

    int norbA_;
    int norbB_;

  public:
    DimerDipole(const std::shared_ptr<const Reference> ref, const int nstart, const int nfenceA, const int nfenceB,
      std::shared_ptr<const Coeff> coeff) : CIDipole(ref, nstart, nfenceB, coeff), norbA_(nfenceA-nstart), norbB_(nfenceB-nfenceA)
    {

      for (int i = 0; i < 3; ++i) {
        dipole_moA_[i] = dipole_mo_[i]->get_submatrix(nocc_, nocc_, norbA_, norbA_);
        dipole_moB_[i] = dipole_mo_[i]->get_submatrix(nocc_ + norbA_, nocc_ + norbA_, norbB_, norbB_);
        cross_dipole_[i] = dipole_mo_[i]->get_submatrix(nocc_, nocc_ + norbA_, norbA_, norbB_);
      }
    }

  template <int unit> std::shared_ptr<const Matrix> dipoles(const int i) const { return (unit==0? dipole_moA_[i] : dipole_moB_[i]); }
  std::shared_ptr<const Matrix> cross_dipole(const int i) const { return cross_dipole_[i]; }

  void compute(std::shared_ptr<const Dvec>) override { assert(false); } // disable
  void print() const override { assert(false); } // disable

};

}

#endif
