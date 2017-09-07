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

#include <src/wfn/hcoreinfo.h>
#include <src/wfn/geometry.h>

// Contains info on DKH2 analytic gradient.
//

namespace bagel {

// metadata
struct DKH2AnalyticData {
    int natom;
    int nbasis;
    int nunc;

    std::shared_ptr<const Molecule> molu;

    Matrix U_T;
    std::vector<Matrix> PU;

    std::vector<Matrix> s_X;
    std::vector<Matrix> T_pX;
    std::vector<Matrix> V_pX;
    std::vector<Matrix> O_pX;

    Matrix id;
    std::map<std::shared_ptr<const VectorB>, std::shared_ptr<Matrix>> vec2mat;
};

class DKH2Analytic : public HcoreInfo {
  private:
    std::shared_ptr<DKH2AnalyticData> gradinit(std::shared_ptr<const Geometry>) const;
    void contracts(std::shared_ptr<const Geometry>, std::shared_ptr<const Matrix>, std::shared_ptr<DKH2AnalyticData>) const;
    void overlapgrad(std::shared_ptr<const Geometry>, std::shared_ptr<const Matrix>, std::shared_ptr<DKH2AnalyticData>) const;
    void kineticgrad(std::shared_ptr<const Geometry>, std::shared_ptr<const Matrix>, std::shared_ptr<DKH2AnalyticData>) const;
    void naigrad(std::shared_ptr<const Geometry>, std::shared_ptr<const Matrix>, std::shared_ptr<DKH2AnalyticData>) const;
    void smallnaigrad(std::shared_ptr<const Geometry>, std::shared_ptr<const Matrix>, std::shared_ptr<DKH2AnalyticData>) const;
    void store_mat(std::shared_ptr<const VectorB>, std::shared_ptr<DKH2AnalyticData>) const;

  public:
    DKH2Analytic() : HcoreInfo() { }
    DKH2Analytic(std::shared_ptr<const PTree> idata) : HcoreInfo(idata) { }

    std::vector<std::shared_ptr<Matrix>> dkh_grad(std::shared_ptr<const Molecule>) const;
};

}

#endif
