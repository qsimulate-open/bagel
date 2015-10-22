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

#include <src/util/math/sphharmonics.h>
#include <src/util/math/gamma.h>
#include <src/integral/rys/erirootlist.h>
#include <src/util/math/factorial.h>
#include <src/periodic/simulationcell.h>
#include <src/periodic/pdata.h>
#include <src/util/parallel/resources.h>
#include <src/periodic/tree.h>

namespace bagel {

class PFMM {
  protected:
    std::shared_ptr<const SimulationCell> scell_;
    const bool dodf_;
    int lmax_, ws_;
    int extent_sum_;
    double thresh_;
    int ndim_;
    int msize_, osize_; // #multipoles in M and O
    std::vector<std::complex<double>> mlm_;

    // Mlm
    int max_rank_;
    double *rvec_, *kvec_;
    double* T_;
    double* Rsq_;
    double *roots_, *weights_;

    // near-field FMM
    int max_height_;
    bool do_contract_;

    double dot(std::array<double, 3> b, std::array<double, 3> c) { return b[0]*c[0]+b[1]*c[1]+b[2]*c[2]; }
    std::array<double, 3> cross(std::array<double, 3> b, std::array<double, 3> c, double s = 1.0) {
      std::array<double, 3> out;
      out[0] = (b[1]*c[2] - b[2]*c[1]) * s;
      out[1] = (b[2]*c[0] - b[0]*c[2]) * s;
      out[2] = (b[0]*c[1] - b[1]*c[0]) * s;
      return out;
    }

    bool allocated_here_;
    std::shared_ptr<StackMem> stack_;
    size_t size_allocated_;
    double* buff_;
    std::vector<std::array<int, 3>> generate_vidx(const int n) const;
    void compute_Mlm();
    void compute_Mlm_direct();
    void root_weight(const int l, const int size);
    void allocate_arrays(const size_t ps);
    std::vector<std::complex<double>> compute_Slm(std::shared_ptr<const PData> density) const;
    std::shared_ptr<const PData> compute_far_field(std::shared_ptr<const PData> density) const;
    std::shared_ptr<const PData> compute_cfmm(std::shared_ptr<const PData> density) const;

  public:
    PFMM(std::shared_ptr<const SimulationCell>, const bool dodf = true, const int lmax = 10, const int ws = 2, const int extent = 10,
         const double thresh = PRIM_SCREEN_THRESH, std::shared_ptr<StackMem> stack = nullptr);
    ~PFMM() { }

    static bool sort_vector(std::array<int, 3> v1, std::array<int, 3> v2) {
      int rad1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
      int rad2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
      return rad1 < rad2;
    };

    int lmax() const { return lmax_; }
    int max_rank() const { return max_rank_; }
    int ws() const { return ws_; }
    int extent_sum() const { return extent_sum_; }
    std::complex<double> mlm(const int i) const { return mlm_[i]; }
    std::vector<std::complex<double>> mlm() const { return mlm_; }

    bool is_in_cff(std::array<double, 3> lvector);
    std::shared_ptr<const PData> pcompute_Jop(std::shared_ptr<const PData> density) const;

};

}

#endif
