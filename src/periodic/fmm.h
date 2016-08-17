//
// BAGEL - Parallel electron correlation program.
// Filename: fmm.h
// Copyright (C) 2016 Toru Shiozaki
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


#ifndef __SRC_PERIODIC_FMM_H
#define __SRC_PERIODIC_FMM_H

#include <set>
#include <src/wfn/geometry.h>
#include <src/periodic/box.h>

namespace bagel {

class Box;
class FMM {
  protected:
    std::shared_ptr<const Geometry>geom_;
    int ns_;
    int lmax_;
    int nbox_, nsp_, nleaf_, nbasis_;
    std::vector<std::array<double, 3>> coordinates_;
    std::array<double, 3> centre_;
    std::array<double, 3> maxxyz_;
    double boxsize_, unitsize_;
    std::vector<int> nbranch_;

    std::vector<std::shared_ptr<Box>> box_;
    double thresh_;
    int ws_;

    void init();
    void get_boxes();

  public:
    FMM(std::shared_ptr<const Geometry> geom, const int ns, const int lmax = 10, const double thresh = PRIM_SCREEN_THRESH, const int ws = 0);
    ~FMM() { }

    std::array<double, 3> centre() const { return centre_; }
    double boxsize() const { return boxsize_; }
    void M2M();
    void M2L();
    void L2L();
    std::shared_ptr<const ZMatrix> compute_energy(std::shared_ptr<const Matrix> density = nullptr, const bool dodf = false, const double schwarz_thresh = 0.0) const;
    void print_boxes(const int i) const;
};

}
#endif
