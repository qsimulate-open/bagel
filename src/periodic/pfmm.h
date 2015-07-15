//
// BAGEL - Parallel electron correlation program.
// Filename: pfmm.h
// Copyright (C) 2015 Toru Shiozaki
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


#ifndef __SRC_PERIODIC_PFMM_H
#define __SRC_PERIODIC_PFMM_H

#include <src/periodic/simulationcell.h>

namespace bagel {

class PFMM {
  protected:
    std::shared_ptr<const SimulationCell> scell_;
    int lmax_, ws_;
    int num_multipoles_;

  public:
    PFMM(std::shared_ptr<const SimulationCell>, const int lmax = ANG_HRR_END, const int ws = 2);
    ~PFMM() { }

    bool is_in_cff(std::array<double, 3> lvector);
    void print_multipoles(std::vector<std::complex<double>> multipoles) const;

};

}

#endif
