//
// BAGEL - Parallel electron correlation program.
// Filename: denom.h
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

#ifndef __SRC_SMITH_DENOM_H
#define __SRC_SMITH_DENOM_H

#include <src/util/matrix.h>
#include <src/wfn/rdm.h>

namespace bagel {
namespace SMITH {

class Denom {
  protected:
    std::shared_ptr<const Matrix> shalf_hh_; 
    std::shared_ptr<const Matrix> shalf_xh_; 
    std::shared_ptr<const Matrix> shalf_xhh_; 
    std::shared_ptr<const Matrix> shalf_xxh_; 

    std::unique_ptr<double[]> denom_hh_;
    std::unique_ptr<double[]> denom_xh_;
    std::unique_ptr<double[]> denom_xhh_;
    std::unique_ptr<double[]> denom_xxh_;

    // init functions
    void init_hh_(const RDM<1>&, const RDM<2>&, const RDM<3>&, const RDM<4>&, const Matrix& fock);
    void init_xh_(const RDM<1>&, const RDM<2>&, const RDM<3>&, const RDM<4>&, const Matrix& fock);
    void init_xhh_(const RDM<1>&, const RDM<2>&, const RDM<3>&, const RDM<4>&, const Matrix& fock);
    void init_xxh_(const RDM<1>&, const RDM<2>&, const RDM<3>&, const RDM<4>&, const Matrix& fock);

  public:
    Denom(const RDM<1>&, const RDM<2>&, const RDM<3>&, const RDM<4>&, const Matrix& fock);

    std::shared_ptr<const Matrix> shalf_hh() const { return shalf_hh_; } 
    std::shared_ptr<const Matrix> shalf_xh() const { return shalf_xh_; } 
    std::shared_ptr<const Matrix> shalf_xhh() const { return shalf_xhh_; } 
    std::shared_ptr<const Matrix> shalf_xxh() const { return shalf_xxh_; } 

    const double& denom_hh(const size_t i) const { return denom_hh_[i]; } 
    const double& denom_xh(const size_t i) const { return denom_xh_[i]; } 
    const double& denom_xhh(const size_t i) const { return denom_xhh_[i]; } 
    const double& denom_xxh(const size_t i) const { return denom_xxh_[i]; } 

};

}
}

#endif
