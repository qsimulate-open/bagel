//
// BAGEL - Parallel electron correlation program.
// Filename: cpcasscf.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_GRAD_CPCASSCF_H
#define __SRC_GRAD_CPCASSCF_H

#include <src/util/matrix.h>
#include <memory>
#include <src/wfn/reference.h>
#include <src/fci/civec.h>
#include <src/fci/fci.h>
#include <src/util/pairfile.h>

// CP-CASSCF Z vector equation. Note that in addition to orbital derivatives (as in CPHF),
// we have CI derivatives in CP-CASSCF.

namespace bagel {

class CPCASSCF {
  protected:
    const std::shared_ptr<const PairFile<Matrix, Dvec>> grad_;
    const std::shared_ptr<const Dvec> civector_;
    const std::shared_ptr<const Matrix> eig_;
    const std::shared_ptr<const DFHalfDist> half_;
    const std::shared_ptr<const DFHalfDist> halfjj_;
    const std::shared_ptr<const Reference> ref_;
    const std::shared_ptr<const Geometry> geom_;
    const std::shared_ptr<const FCI> fci_;


    std::shared_ptr<Matrix> compute_amat(std::shared_ptr<const Dvec> z1, std::shared_ptr<const Dvec> c1, std::shared_ptr<const Determinants>) const;

  public:
    CPCASSCF(const std::shared_ptr<const PairFile<Matrix, Dvec>> grad, const std::shared_ptr<const Dvec> c, const std::shared_ptr<const Matrix> eig,
             const std::shared_ptr<const DFHalfDist> half, const std::shared_ptr<const DFHalfDist> halfjj,
             const std::shared_ptr<const Reference> g, const std::shared_ptr<const FCI> f);

    std::shared_ptr<PairFile<Matrix, Dvec>> solve() const;

};

}

#endif

