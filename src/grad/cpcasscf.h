//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cpcasscf.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_GRAD_CPCASSCF_H
#define __SRC_GRAD_CPCASSCF_H

#include <src/ci/fci/distfci.h>
#include <src/ci/fci/fci.h>
#include <src/ci/fci/knowles.h>
#include <src/util/math/pairfile.h>
#include <src/multi/casscf/qvec.h>

// CP-CASSCF Z vector equation. Note that in addition to orbital derivatives (as in CPHF),
// we have CI derivatives in CP-CASSCF.

namespace bagel {

class CPCASSCF {
  protected:
    std::shared_ptr<const PairFile<Matrix, Dvec>> grad_;
    std::shared_ptr<const Dvec> civector_;

    // half transformed DF integrals (closed+act, J^-1/2 multiplied)
    std::shared_ptr<const DFHalfDist> halfj_;
    std::shared_ptr<const Qvec> qvec_;
    std::shared_ptr<const Matrix> fock_;
    std::shared_ptr<const Matrix> fockinact_;

    std::shared_ptr<const Reference> ref_;
    std::shared_ptr<const Geometry> geom_;

    std::shared_ptr<FCI_base> fci_native_;
    std::shared_ptr<FCI> fci_;

    int ncore_;
    bool imag_;
    std::shared_ptr<const Matrix> coeff_;

    std::shared_ptr<PairFile<Matrix,Dvec>> form_sigma(std::shared_ptr<const PairFile<Matrix,Dvec>> z, std::shared_ptr<const DFFullDist>,
                                                      std::shared_ptr<const Determinants> det, std::shared_ptr<const Matrix>,
                                                      const bool antisym, const double lambda) const;

    std::shared_ptr<Matrix> compute_amat(std::shared_ptr<const Dvec> z1, std::shared_ptr<const Dvec> c1, std::shared_ptr<const Determinants>) const;
    std::tuple<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>,std::shared_ptr<Matrix>> compute_orb_denom_and_fock() const;

  public:
    CPCASSCF(std::shared_ptr<const PairFile<Matrix, Dvec>> grad, std::shared_ptr<const Dvec> c, std::shared_ptr<const DFHalfDist> halfj,
             std::shared_ptr<const Reference> g, std::shared_ptr<FCI_base> f, const int ncore = 0, const bool imag = false, std::shared_ptr<const Matrix> coeff = nullptr);

    // tuple of Z, z, and X.
    std::tuple<std::shared_ptr<const Matrix>, std::shared_ptr<const Dvec>, std::shared_ptr<const Matrix>, std::shared_ptr<const Matrix>>
      solve(const double thresh, const int maxiter = 100, std::shared_ptr<const Matrix> den = nullptr, const bool project_all = false);

};

}

#endif

