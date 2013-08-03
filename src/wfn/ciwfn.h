//
// BAGEL - Parallel electron correlation program.
// Filename: ciwfn.h
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


#ifndef _BAGEL_WFN_CIWFN_H
#define _BAGEL_WFN_CIWFN_H

#include <src/scf/coeff.h>
#include <src/wfn/geometry.h>
#include <src/fci/dvec.h>

// Stores the result of some CI type wavefunction (FCI, CASSCF, etc.)
// Right now only contains the bare minimum needed in the Dimer class

namespace bagel {

class CIWfn {

  protected:
    // Geometry which this wave function is belonging to
    const std::shared_ptr<const Geometry> geom_;
    // MO coefficients
    const std::shared_ptr<const Coeff> coeff_;

    const int ncore_;
    const int nact_;
    const int nvirt_;

    const int nstates_;
    std::shared_ptr<const Determinants> det_;
    std::shared_ptr<const Dvec> ccvec_;
    std::vector<double> energies_;

  public:
    CIWfn(std::shared_ptr<const Geometry> g, std::shared_ptr<const Coeff> c,
              const int ncore, const int nact, const int nvirt,
              std::vector<double> en, std::shared_ptr<const Dvec> ccvec)
      : geom_(g), coeff_(c), ncore_(ncore), nact_(nact), nvirt_(nvirt), nstates_(ccvec->ij()),
         det_(ccvec->det()), ccvec_(ccvec), energies_(en) {}

    std::shared_ptr<const Geometry> geom() const { return geom_; }
    const std::shared_ptr<const Coeff> coeff() const { return coeff_; }

    int ncore() const { return ncore_; }
    int nact() const { return nact_; }
    int nvirt() const {return nvirt_; }

    int nstates() const { return nstates_; }

    std::vector<double> energies() const { return energies_; }
    double energy(int i) const {return energies_[i];}

    // function to return a CI vectors from orbital info
    std::shared_ptr<const Determinants> det() const { return det_; }
    std::shared_ptr<const Dvec> civectors() const { return ccvec_; }
};

}

#endif
