//
// BAGEL - Parallel electron correlation program.
// Filename: dfdist.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#ifndef __SRC_DF_DFDIST_H
#define __SRC_DF_DFDIST_H

#include <src/df/df.h>
#include <src/rysint/eribatch.h>
#include <src/rysint/libint.h>

namespace bagel {

class DFHalfDist;

class DFDist : public DensityFit {
  protected:
    std::list<std::shared_ptr<DFBlock> > blocks_;

    std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) override;

    void common_init(const std::vector<std::shared_ptr<const Atom> >&,
                     const std::vector<std::shared_ptr<const Atom> >&,
                     const std::vector<std::shared_ptr<const Atom> >&, const double thresh, const bool compute_inv) override;

  public:
    // construction of a block from AO integrals
    DFDist(const int nbas, const int naux, const std::vector<std::shared_ptr<const Atom> >& atoms,
                                           const std::vector<std::shared_ptr<const Atom> >& aux_atoms, const double thr, const bool inverse)
      : DensityFit(nbas, naux) {
      common_init(atoms, atoms, aux_atoms, thr, inverse);
    }

    // compute half transforms; c is dimensioned by nbasis_;
    std::shared_ptr<DFHalfDist> compute_half_transform(const double* c, const size_t nocc) const;

    // compute a J operator, given density matrices in AO basis
    std::unique_ptr<double[]> compute_Jop(const double* den) const override;

    std::unique_ptr<double[]> compute_cd(const double* den) const override;

};


class DFHalfDist {

};

}


#endif
