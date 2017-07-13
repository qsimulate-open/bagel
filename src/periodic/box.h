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


#ifndef __SRC_PERIODIC_BOX_H
#define __SRC_PERIODIC_BOX_H

#include <src/util/constants.h>
#include <src/util/parallel/resources.h>
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
    std::shared_ptr<const Box> parent_;
    std::vector<std::weak_ptr<const Box>> child_;
    std::vector<std::shared_ptr<const Box>> inter_;
    std::vector<std::shared_ptr<const Box>> neigh_;
    std::vector<std::shared_ptr<const Box>> nonneigh_;
    std::vector<std::shared_ptr<const ShellPair>> sp_;

    double extent_, schwarz_thresh_;
    int nbasis0_, nbasis1_;
    int nmult_;
    int nocc_;
    int nsp_;
    size_t nsize_, msize_, olm_ndim_, olm_mdim_, olm_size_block_;

    std::vector<std::shared_ptr<ZMatrix>> box_olm_;
    std::shared_ptr<ZMatrix> olm_ji_;
    std::shared_ptr<ZMatrix> mlm_ji_;
    std::vector<std::array<int, 2>> offsets_;
    std::vector<std::pair<int, int>> coffsets_s_, coffsets_u_;

    void init();
    void insert_sp(const std::vector<std::shared_ptr<const ShellPair>>&);
    void insert_child(std::shared_ptr<const Box> = NULL);
    void insert_parent(std::shared_ptr<const Box> parent = NULL);
    bool is_neigh(std::shared_ptr<const Box> b, const double ws) const;
    void get_neigh(const std::vector<std::shared_ptr<Box>>& box, const double ws);
    void get_inter(const std::vector<std::shared_ptr<Box>>& box, const double ws);

    std::shared_ptr<ZVectorB> multipole_;
    std::shared_ptr<ZVectorB> localJ_;
    void compute_M2M(std::shared_ptr<const Matrix> density);
    void compute_M2M_X(std::shared_ptr<const Matrix> ocoeff_sj, std::shared_ptr<const Matrix> ocoeff_ui);
    void compute_multipolesX();
    void sort_sp();
    std::shared_ptr<const ZVectorB> shift_multipoles(std::shared_ptr<const ZVectorB> oa, std::array<double, 3> rab) const;
    std::shared_ptr<const ZMatrix> shift_multipolesX(const int lmax, std::shared_ptr<const ZMatrix> oa, std::array<double, 3> rab) const;
    std::shared_ptr<const ZVectorB> shift_localL(std::shared_ptr<const ZVectorB> mr, std::array<double, 3> rb) const;
    std::shared_ptr<const ZMatrix> shift_localLX(const int lmax, std::shared_ptr<const ZMatrix> mr, std::array<double, 3> rb) const;
    std::shared_ptr<const ZVectorB> shift_localM(std::shared_ptr<const ZVectorB> olm, std::array<double, 3> r12) const;
    std::shared_ptr<const ZMatrix> shift_localMX(const int lmax, std::shared_ptr<const ZMatrix> olm, std::array<double, 3> r12) const;
    void compute_M2L();
    void compute_M2L_X();
    void compute_L2L();
    void compute_L2L_X();
    std::shared_ptr<const Matrix> compute_exact_ff(std::shared_ptr<const Matrix> density) const; //debug
    std::shared_ptr<const Matrix> compute_exact_ff_K(std::shared_ptr<const Matrix> density) const; //debug
    std::shared_ptr<const Matrix> compute_exact_ff_J(std::shared_ptr<const Matrix> density) const; //debug
    std::shared_ptr<const Matrix> compute_Fock_nf(std::shared_ptr<const Matrix> density, std::shared_ptr<const VectorB> max_den) const;
    std::shared_ptr<const Matrix> compute_Fock_ff(std::shared_ptr<const Matrix> density) const;
    std::shared_ptr<const Matrix> compute_Fock_ff_K(std::shared_ptr<const Matrix> ocoeff_ti) const;
    // temporary: allow constructing FMM_J and FMM_K separately with different parameters
    std::shared_ptr<const Matrix> compute_Fock_nf_J(std::shared_ptr<const Matrix> density, std::shared_ptr<const VectorB> max_den) const;
    std::shared_ptr<const Matrix> compute_Fock_nf_K(std::shared_ptr<const Matrix> density, std::shared_ptr<const VectorB> max_den) const;
    double coulomb_ff() const;
    double exchange_ff() const;

  public:
    Box(int n, double size, const std::array<double, 3>& c, const int id, const std::array<int, 3>& v, const int lmax = 10,
        const int lmax_k = 10, const std::vector<std::shared_ptr<const ShellPair>>& sp = std::vector<std::shared_ptr<const ShellPair>>(),
        const double schwarz = 0.0)
     : rank_(n), boxsize_(size), centre_(c), boxid_(id), tvec_(v), lmax_(lmax), lmax_k_(lmax_k), sp_(sp), schwarz_thresh_(schwarz) { }

    Box(std::shared_ptr<const Box> b, const std::vector<std::shared_ptr<const ShellPair>>& sp)
     : rank_(b->rank()), boxsize_(b->boxsize()), centre_(b->centre()), boxid_(b->boxid()), tvec_(b->tvec()), lmax_(b->lmax()), lmax_k_(b->lmax_k()),
       sp_(sp), schwarz_thresh_(b->schwarz_thresh()) { }

    ~Box() { }

    const std::array<double, 3>& centre() const { return centre_; }
    double centre(const int i) const { return centre_[i]; }
    double boxsize() const { return boxsize_; }

    int rank() const { return rank_; }
    int boxid() const { return boxid_; }
    const std::array<int, 3>& tvec() const { return tvec_; }
    int lmax() const { return lmax_; }
    int lmax_k() const { return lmax_k_; }
    double schwarz_thresh() const { return schwarz_thresh_; }

    int nsp() const { return nsp_; }
    int nchild() const { return child_.size(); }
    int nneigh() const { return neigh_.size(); }
    int ninter() const { return inter_.size(); }
    double extent() const { return extent_; }
    int nbasis0() const { return nbasis0_; }
    int nbasis1() const { return nbasis1_; }

    std::shared_ptr<const Box> parent() const { return parent_; }
    std::shared_ptr<const Box> child(const int i) const { return child_[i].lock(); }
    const std::vector<std::shared_ptr<const Box>>& neigh() const { return neigh_; }
    const std::vector<std::shared_ptr<const Box>>& interaction_list() const { return inter_; }
    const std::vector<std::shared_ptr<const Box>>& nonneigh() const { return nonneigh_; }
    const std::vector<std::shared_ptr<const ShellPair>>& sp() const { return sp_; }
    std::shared_ptr<const ShellPair> sp(const int i) const { return sp_[i]; }
    const std::array<std::shared_ptr<const Shell>, 2>& shells(const int i) const { return sp_[i]->shells(); }

    std::shared_ptr<ZVectorB> multipole() const { return multipole_; }
    std::shared_ptr<ZVectorB> localJ() const { return localJ_; }

    std::shared_ptr<const ZMatrix> olm_ji() const { return olm_ji_; }
    std::shared_ptr<const ZMatrix> mlm_ji() const { return mlm_ji_; }
    const std::vector<std::shared_ptr<ZMatrix>>& box_olm() const { return box_olm_; }
    int box_olm_size() const { return box_olm_[0]->size(); }
    std::vector<std::shared_ptr<ZMatrix>>& box_olm() { return box_olm_; }
    std::shared_ptr<ZMatrix> box_olm(const int i) const { return box_olm_[i]; }

    void print_box() const;
};

}
#endif
