//
// BAGEL - Parallel electron correlation program.
// Filename: ecpbatch.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_INTEGRAL_ECP_ECPBATCH_H
#define __SRC_INTEGRAL_ECP_ECPBATCH_H

#include <tuple>
#include <src/molecule/molecule.h>
#include <src/util/constants.h>
#include <src/integral/integral.h>

namespace bagel {

class ECPBatch: public Integral_base<double> {
  protected:

    std::shared_ptr<const Shell_ECP> ecp_basisinfo_;
    std::array<std::shared_ptr<const Shell>,2> basisinfo_;
    std::shared_ptr<const Molecule> mol_;


  public:
    ECPBatch(const std::shared_ptr<const Shell_ECP>& ecp_info, const std::array<std::shared_ptr<const Shell>,2>& info,
            const std::shared_ptr<const Molecule> mol)
     : ecp_basisinfo_(ecp_info), basisinfo_(info), mol_(mol) {
//   const double integral_thresh = PRIM_SCREEN_THRESH;
      /// call some functions to compute integrals
    }

    ~ECPBatch() {}

    void compute() override;

};

}

#endif
