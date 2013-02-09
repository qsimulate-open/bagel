//
// BAGEL - Parallel electron correlation program.
// Filename: dfdata.h
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


#ifndef __SRC_REL_DFDATA_H
#define __SRC_REL_DFDATA_H

#include <memory>
#include <string>
#include <map>
#include <src/df/df.h>
#include <src/wfn/reference.h>
#include <src/rel/alpha.h>
#include <src/util/zmatrix.h>

namespace bagel {

class DFData {
  protected:
    std::shared_ptr<const DFDist> dfdata_;
    std::pair<int, int> coord_;
    std::pair<int, int> basis_; 
    bool swap_;
    std::shared_ptr<Alpha> alpha_;
    std::shared_ptr<Sigma> sigma1_;
    std::shared_ptr<Sigma> sigma2_;
    std::pair<std::shared_ptr<ZMatrix>, std::shared_ptr<ZMatrix>> spinor_;
    std::complex<double> fac_;

    std::pair<std::complex<double>, std::complex<double>> calc_coeff(std::pair<const int, const int>, std::pair<const int, const int>); 

    DFData(const DFData&, bool , bool);
    std::pair<std::shared_ptr<ZMatrix>, std::shared_ptr<ZMatrix>> compute_spinor(std::pair<const int, const int>, std::pair<const int, const int>);

  public:
    DFData(std::shared_ptr<const DFDist>, std::pair<int, int>, const int);
    DFData(const DFData&) = delete;
    DFData() = delete;

    std::shared_ptr<const DFDist> df() const { return dfdata_; }
    std::pair<int, int> coord() const { return coord_; }
    std::pair<int, int> basis() const { return basis_; }
    bool cross() const { return coord_.first != coord_.second; }
    std::shared_ptr<Alpha> alpha() const { return alpha_; }
    std::complex<double> fac() const { return fac_; }
    bool swapped() const { return swap_; }
    int coeff_index() const;
    std::shared_ptr<const DFData> opp();
    std::shared_ptr<const DFData> swap();
    std::shared_ptr<const DFData> opp_and_swap();

    const std::tuple<int, int, int, int> compute_index_Jop() const;

};

}

#endif
