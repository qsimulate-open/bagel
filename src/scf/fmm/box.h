//
// BAGEL - Parallel electron correlation program.
// Filename: box.h
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


#ifndef __SRC_SCF_FMM_BOX_H
#define __SRC_SCF_FMM_BOX_H

#include <src/molecule/shellpair.h>

namespace bagel {

class Box {
  friend class FMM;
  protected:
    int rank_;
    double boxsize_;
    std::array<double, 3> centre_;
    int boxid_;
    std::array<int, 3> tvec_;
    int lmax_, lmax_k_, ws_;
    std::weak_ptr<const Box> parent_;
    std::vector<std::weak_ptr<const Box>> child_;
    std::vector<std::weak_ptr<const Box>> inter_;
    std::vector<std::weak_ptr<const Box>> neigh_;
    std::vector<std::weak_ptr<const Box>> nonneigh_; //for debugging
    std::vector<std::shared_ptr<const ShellPair>> sp_;
    std::map<std::shared_ptr<const Shell>, std::tuple<int/*local_offset*/,int/*all_offset*/,int/*num*/>> shell0_;
    int nchild_, ninter_, nneigh_;

    double extent_, schwarz_thresh_;
    int nmult_;
    int nsp_;
    int nshell0_;
    size_t nsize_, msize_, olm_ndim_, olm_mdim_, olm_size_block_;

    std::shared_ptr<ZMatrix> olm_ji_;
    std::shared_ptr<ZMatrix> mlm_ji_;

    void init();
    void insert_sp(const std::vector<std::shared_ptr<const ShellPair>>&);
    void insert_child(std::weak_ptr<const Box>);
    void insert_parent(std::weak_ptr<const Box> parent);
    bool is_neigh(std::weak_ptr<const Box> b, const double ws) const;
    void get_neigh(const std::vector<std::weak_ptr<Box>>& box, const double ws);
    void get_inter(const std::vector<std::weak_ptr<Box>>& box, const double ws);
    void sort_sp();

    std::shared_ptr<ZVectorB> olm_;
    std::shared_ptr<ZVectorB> mlm_;

    void compute_M2M(std::shared_ptr<const Matrix> density);
    void compute_M2M_X(std::shared_ptr<const Matrix> ocoeff_sj, std::shared_ptr<const Matrix> ocoeff_ui);
    void compute_M2L();
    void compute_M2L_X();
    void compute_L2L();
    void compute_L2L_X();
    std::shared_ptr<const ZMatrix> shift_multipolesX(const int lmax, std::shared_ptr<const ZMatrix> oa, std::array<double, 3> rab) const;
    std::shared_ptr<const ZMatrix> shift_localLX(const int lmax, std::shared_ptr<const ZMatrix> mr, std::array<double, 3> rb) const;
    std::shared_ptr<const ZMatrix> shift_localMX(const int lmax, std::shared_ptr<const ZMatrix> olm, std::array<double, 3> r12) const;

    std::shared_ptr<const Matrix> compute_exact_ff(std::shared_ptr<const Matrix> density) const; //debug
    std::shared_ptr<const Matrix> compute_Fock_nf(std::shared_ptr<const Matrix> density, std::shared_ptr<const VectorB> max_den) const;
    std::shared_ptr<const Matrix> compute_Fock_ff(std::shared_ptr<const Matrix> density) const;
    std::shared_ptr<const Matrix> compute_Fock_ff_K(std::shared_ptr<const Matrix> ocoeff_ti) const;
    // allow constructing FMM_J and FMM_K separately with different parameters
    std::shared_ptr<const Matrix> compute_Fock_nf_J(std::shared_ptr<const Matrix> density, std::shared_ptr<const VectorB> max_den) const;
    std::shared_ptr<const Matrix> compute_Fock_nf_K(std::shared_ptr<const Matrix> density, std::shared_ptr<const VectorB> max_den) const;

  public:
    Box() { }
    Box(int n, double size, const std::array<double, 3>& c, const int id, const std::array<int, 3>& v, const int lmax = 10,
        const int lmax_k = 10, const std::vector<std::shared_ptr<const ShellPair>>& sp = std::vector<std::shared_ptr<const ShellPair>>(),
        const double schwarz = 0.0)
     : rank_(n), boxsize_(size), centre_(c), boxid_(id), tvec_(v), lmax_(lmax), lmax_k_(lmax_k), sp_(sp), schwarz_thresh_(schwarz) { }

    ~Box() { }

    const std::array<double, 3>& centre() const { return centre_; }
    double centre(const int i) const { return centre_[i]; }
    double boxsize() const { return boxsize_; }

    int rank() const { return rank_; }
    int boxid() const { return boxid_; }
    const std::array<int, 3>& tvec() const { return tvec_; }
    double extent() const { return extent_; }
    std::shared_ptr<const Box> parent() const { return parent_.lock(); }

    const std::vector<std::shared_ptr<const ShellPair>>& sp() const { return sp_; }
    const std::vector<std::weak_ptr<const Box>>& child() const { return child_; }

    std::shared_ptr<ZVectorB> olm() const { return olm_; }
    std::shared_ptr<ZVectorB> mlm() const { return mlm_; }

    std::shared_ptr<const ZMatrix> olm_ji() const { return olm_ji_; }
    std::shared_ptr<const ZMatrix> mlm_ji() const { return mlm_ji_; }

    void print_box() const;
};

}
#endif
