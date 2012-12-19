//
// BAGEL - Parallel electron correlation program.
// Filename: eribatch_base.h
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


#ifndef __SRC_RYSINT_ERIBATCH_BASE_H
#define __SRC_RYSINT_ERIBATCH_BASE_H

#include <src/rysint/macros.h>
#include <src/scf/shell.h>
#include <src/rysint/rysint.h>
#include <tuple>

namespace bagel {

class ERIBatch_base : public RysInt{
  protected:
    const bool breit_;

    void root_weight(const int ps);
    void compute_ssss(const double);

  public:
    ERIBatch_base(const std::array<std::shared_ptr<const Shell>,4>& o, const double max_density, const int deriv, const bool breit = false) : RysInt(o), breit_(breit) {
      const double integral_thresh = (max_density != 0.0) ? (PRIM_SCREEN_THRESH / max_density) : 0.0;
      deriv_rank_ = deriv;

      // determins if we want to swap shells
      set_swap_info(true);

      // stores AB and CD
      set_ab_cd();

      // set primsize_ and contsize_, as well as relevant members
      set_prim_contsizes();

      // sets angular info
      int asize_final, csize_final, asize_final_sph, csize_final_sph;
      std::tie(asize_final, csize_final, asize_final_sph, csize_final_sph) = set_angular_info();

      // allocate
      allocate_data(asize_final, csize_final, asize_final_sph, csize_final_sph);
      allocate_arrays(primsize_);

      compute_ssss(integral_thresh);

      root_weight(primsize_);

    };
    ~ERIBatch_base() {};

};

}

#endif
