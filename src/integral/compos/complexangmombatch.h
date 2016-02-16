//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexangmombatch.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

// Orbital angular momentum integrals, for GIAO basis functions, around some center mcoord_
// We return <A| L |B>; to agree with the zero-field code (AngMomBatch) one or the other has to be scaled

#ifndef __SRC_INTEGRAL_COMPOS_COMPLEXANGMOMBATCH_H
#define __SRC_INTEGRAL_COMPOS_COMPLEXANGMOMBATCH_H

#include <src/integral/os/osintegral.h>
#include <src/molecule/molecule.h>

namespace bagel {

class ComplexAngMomBatch : public OSIntegral<std::complex<double>,Int_t::London> {
  protected:
    std::array<double,3> magnetic_field_;
    std::array<double,3> mcoord_;

    void perform_VRR(std::complex<double>*) override;
    virtual std::complex<double> get_P(const double coord1, const double coord2, const double exp1, const double exp2, const double one12,
                                       const int dim, const bool swap) override;

    int nblocks() const override { return 3; }
    int nrank() const override { return 0; }

  public:
    ComplexAngMomBatch(const std::array<std::shared_ptr<const Shell>,2>& basis, const std::array<double,3> _magnetic_field, const std::array<double,3> _mcoord)
      : OSIntegral<std::complex<double>,Int_t::London>(basis), magnetic_field_(_magnetic_field), mcoord_(_mcoord) { common_init(); }

    constexpr static int Nblocks() { return 3; }

    void compute() override;
};

}

#endif
