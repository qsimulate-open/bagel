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


#ifndef __SRC_GRAD_CPCASSCF_H
#define __SRC_GRAD_CPCASSCF_H

#include <src/fci/fci.h>
#include <src/math/pairfile.h>

// CP-CASSCF Z vector equation. Note that in addition to orbital derivatives (as in CPHF),
// we have CI derivatives in CP-CASSCF.

namespace bagel {

class CPCASSCF {
  protected:
    std::shared_ptr<const PairFile<Matrix, Dvec>> grad_;
    std::shared_ptr<const Dvec> civector_;
    std::shared_ptr<const DFHalfDist> half_;
    std::shared_ptr<const DFHalfDist> halfjj_;
    std::shared_ptr<const Reference> ref_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const FCI> fci_;
    std::shared_ptr<const Matrix> coeff_;

    std::shared_ptr<PairFile<Matrix,Dvec>> form_sigma(std::shared_ptr<const PairFile<Matrix,Dvec>> z, std::shared_ptr<const DFHalfDist>,
                                                      std::shared_ptr<const DFFullDist>, std::shared_ptr<const Determinants> det, std::shared_ptr<const Matrix>) const;
    std::shared_ptr<Matrix> form_sigma_sym(std::shared_ptr<const PairFile<Matrix,Dvec>> z, std::shared_ptr<const DFHalfDist>,
                                           std::shared_ptr<const DFFullDist>, std::shared_ptr<const Determinants> det, std::shared_ptr<const Matrix>) const;
    std::shared_ptr<Matrix> compute_amat(std::shared_ptr<const Dvec> z1, std::shared_ptr<const Dvec> c1, std::shared_ptr<const Determinants>) const;
    std::shared_ptr<Matrix> compute_orb_denom() const;

  public:
    CPCASSCF(std::shared_ptr<const PairFile<Matrix, Dvec>> grad, std::shared_ptr<const Dvec> c,
             std::shared_ptr<const DFHalfDist> half, std::shared_ptr<const DFHalfDist> halfjj,
             std::shared_ptr<const Reference> g, std::shared_ptr<const FCI> f, std::shared_ptr<const Matrix> coeff = nullptr);

    // tuple of Z, z, and X.
    std::tuple<std::shared_ptr<const Matrix>, std::shared_ptr<const Dvec>, std::shared_ptr<const Matrix>>
      solve(const double thresh, const int maxiter = 100);

};

}

#endif

