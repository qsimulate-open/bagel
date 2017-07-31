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
    std::vector<int> isp_;
    bool do_ff_;
    std::vector<std::array<double, 3>> coordinates_;
    std::array<double, 3> centre_;
    std::array<double, 3> maxxyz_;
    double boxsize_, unitsize_;
    std::vector<int> nbranch_;

    std::vector<std::shared_ptr<Box>> box_;
    double ws_;
    double base_extent_;
    bool do_exchange_, do_cadfj_;
    int lmax_k_;
    bool debug_;
    int xbatchsize_;

    void init();
    void get_boxes();
    void M2M(std::shared_ptr<const Matrix> mat, const bool do_exchange = false) const;
    void M2M_X(std::shared_ptr<const Matrix> ocoeff_sj, std::shared_ptr<const Matrix> ocoeff_ui) const;
    void M2L(const bool do_exchange = false) const;
    void L2L(const bool do_exchange = false) const;

    //NF-CADF-J
    std::shared_ptr<const Matrix> XY_;
    void compute_2index();

  public:
    FMM(std::shared_ptr<const Geometry> geom, const int ns = 4, const int lmax = 10, const double ws = 0.0,
        const bool do_exchange = true, const int lmax_k = 2, const bool do_cadfj = false, const bool debug = true, const int batchsize = -1);
    ~FMM() { }

    const std::array<double, 3>& centre() const { return centre_; }
    double boxsize() const { return boxsize_; }
    bool do_exchange() const { return do_exchange_; }
    std::shared_ptr<const Matrix> compute_Fock_FMM(std::shared_ptr<const Matrix> density = nullptr) const;
    std::shared_ptr<const Matrix> compute_Fock_FMM_K(std::shared_ptr<const Matrix> density = nullptr) const;
    std::shared_ptr<const Matrix> compute_Fock_FMM_J(std::shared_ptr<const Matrix> density = nullptr) const;
    std::shared_ptr<const Matrix> compute_K_ff(std::shared_ptr<const Matrix> ocoeff, std::shared_ptr<const Matrix> overlap) const;
    std::shared_ptr<const Matrix> compute_K_ff_from_den(std::shared_ptr<const Matrix> density, std::shared_ptr<const Matrix> overlap) const;
    void print_boxes(const int i) const;
};

}
#endif
