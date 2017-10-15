//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkh2analytic.h
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#ifndef __SRC_WFN_DKH2ANALYTIC_H
#define __SRC_WFN_DKH2ANALYTIC_H

#include <src/molecule/molecule.h>
#include <src/wfn/diagvec.h>
#include <src/wfn/hcoreinfo.h>

// Contains info on DKH2 analytic gradient.
//

namespace bagel {

class DKH2Analytic : public HcoreInfo {
  private:
    int natom;
    int nunc;

    Matrix U;
    DiagVec s;

    std::vector<Matrix> PU;
    std::vector<Matrix> O_pX;

    Matrix id;

    void gradinit(std::shared_ptr<const Molecule>);
    void smallnaigrad(std::shared_ptr<const Molecule>, std::shared_ptr<const Matrix>);

  public:
    DKH2Analytic() : HcoreInfo() { }
    DKH2Analytic(std::shared_ptr<const PTree> idata) : HcoreInfo(idata) { gradtype_ = false; }

    std::vector<std::shared_ptr<Matrix>> dkh_grad(std::shared_ptr<const Molecule>);
};

}

#endif
