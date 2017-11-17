//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cis.cc
// Copyright (C) 2017 Toru Shiozaki
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

#ifndef __SRC_RESPONSE_CIS_H
#define __SRC_RESPONSE_CIS_H

#include <src/wfn/method.h>

namespace bagel {

// perform CI singles
class CIS : public Method {
  protected:
    const int nstate_;
    const int nocc_;
    const int nvirt_;
    const int maxiter_;

    const double thresh_;

    // orbital energies
    VectorB eig_;
    // CIS energies
    std::vector<double> energy_;

    std::shared_ptr<const DFHalfDist> half_;
    std::shared_ptr<const DFFullDist> fulljj_;
    std::shared_ptr<const Matrix> coeff_; // coeff internally used

    std::vector<std::shared_ptr<const Matrix>> amp_;

  public:
    CIS(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    void compute();
    std::shared_ptr<const Reference> conv_to_ref() const { return ref_; }

    double energy() const { return energy_[0] + ref_->energy(0); }
    std::vector<double> excitation_energy() const { return energy_; }
};

}

#endif
