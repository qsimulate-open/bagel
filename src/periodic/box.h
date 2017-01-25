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
#include <src/molecule/shellpair.h>
#include <src/util/parallel/resources.h>

namespace bagel {

class Box {
  friend class FMM;
  protected:
    int rank_, boxid_, lmax_;
    std::array<int, 3> tvec_;
    std::array<double, 3> centre_;
    double boxsize_;
    std::shared_ptr<const Box> parent_;
    std::vector<std::weak_ptr<const Box>> child_;
    std::vector<std::shared_ptr<const Box>> inter_;
    std::vector<std::shared_ptr<const Box>> neigh_;
    std::vector<std::shared_ptr<const ShellPair>> sp_;

    double thresh_, extent_, schwarz_thresh_;
    int nbasis0_, nbasis1_;
    int nmult_;
    int nocc_;
    int nsize_, msize_;

    std::vector<std::shared_ptr<ZMatrix>> box_olm_;
    std::vector<std::shared_ptr<ZMatrix>> olm_ji_;
    std::vector<std::shared_ptr<ZMatrix>> mlm_ji_;
    std::vector<std::array<int, 2>> offsets_;
    std::vector<std::pair<int, int>> coffsets_s_, coffsets_u_;
    void get_offsets();

    void init();
    void insert_sp(const std::vector<std::shared_ptr<const ShellPair>>&);
    void insert_child(std::shared_ptr<const Box> = NULL);
    void insert_parent(std::shared_ptr<const Box> parent = NULL);
    bool is_neigh(std::shared_ptr<const Box> b, const int ws) const;
    void get_neigh(const std::vector<std::shared_ptr<Box>>& box, const int ws);
    void get_inter(const std::vector<std::shared_ptr<Box>>& box, const int ws);

    std::shared_ptr<ZVectorB> multipole_;
    std::shared_ptr<ZVectorB> localJ_;
    void compute_M2M(std::shared_ptr<const Matrix> density);
    void compute_M2M_X(std::shared_ptr<const Matrix> ocoeff_sj, std::shared_ptr<const Matrix> ocoeff_ui);
    void compute_multipolesX();
    void sort_sp();
    std::shared_ptr<const ZVectorB> shift_multipoles(std::shared_ptr<const ZVectorB> oa, std::array<double, 3> rab) const;
    std::vector<std::shared_ptr<const ZMatrix>> shift_multipolesX(const std::vector<std::shared_ptr<ZMatrix>>& oa, std::array<double, 3> rab) const;
    std::shared_ptr<const ZVectorB> shift_localL(std::shared_ptr<const ZVectorB> mr, std::array<double, 3> rb) const;
    std::vector<std::shared_ptr<const ZMatrix>> shift_localLX(const std::vector<std::shared_ptr<ZMatrix>>& mr, std::array<double, 3> rb) const;
    std::shared_ptr<const ZVectorB> shift_localM(std::shared_ptr<const ZVectorB> olm, std::array<double, 3> r12) const;
    std::vector<std::shared_ptr<const ZMatrix>> shift_localMX(const std::vector<std::shared_ptr<ZMatrix>>& olm, std::array<double, 3> r12) const;
    void compute_M2L();
    void compute_M2L_X();
    void compute_L2L();
    void compute_L2L_X();
    double compute_exact_energy_ff(std::shared_ptr<const Matrix> density) const; //debug
    std::shared_ptr<const ZMatrix> compute_Fock_nf(std::shared_ptr<const Matrix> density, std::vector<double>& max_den) const;
    std::shared_ptr<const ZMatrix> compute_Fock_ff(std::shared_ptr<const Matrix> density) const;
    std::shared_ptr<const ZMatrix> compute_Fock_ffX(std::shared_ptr<const Matrix> ocoeff_ti) const;


  public:
    Box(int n, int id, const std::array<int, 3>& v, const int lmax = 10, const std::vector<std::shared_ptr<const ShellPair>>& sp = std::vector<std::shared_ptr<const ShellPair>>(),
        const double thresh = 0.0)
     : rank_(n), boxid_(id), lmax_(lmax), tvec_(v), sp_(sp), schwarz_thresh_(thresh) { }

    ~Box() { }

    const std::array<double, 3>& centre() const { return centre_; }
    double centre(const int i) const { return centre_[i]; }

    int rank() const { return rank_; }
    int boxid() const { return boxid_; }
    const std::array<int, 3>& tvec() const { return tvec_; }

    int nsp() const { return sp_.size(); }
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
    const std::vector<std::shared_ptr<const ShellPair>>& sp() const { return sp_; }
    std::shared_ptr<const ShellPair> sp(const int i) const { return sp_[i]; }
    const std::array<std::shared_ptr<const Shell>, 2>& shells(const int i) const { return sp_[i]->shells(); }

    std::shared_ptr<ZVectorB> multipole() const { return multipole_; }
    std::shared_ptr<ZVectorB> localJ() const { return localJ_; }

    const std::vector<std::shared_ptr<ZMatrix>>& olm_ji() const { return olm_ji_; }
    std::vector<std::shared_ptr<ZMatrix>>& olm_ji() { return olm_ji_; }
    int olm_ji_size() const { return olm_ji_[0]->size(); }
    const std::vector<std::shared_ptr<ZMatrix>>& mlm_ji() const { return mlm_ji_; }
    std::shared_ptr<ZMatrix> olm_ji(const int i) const { return olm_ji_[i]; }
    std::vector<std::shared_ptr<ZMatrix>>& mlm_ji() { return mlm_ji_; }
    const std::vector<std::shared_ptr<ZMatrix>>& box_olm() const { return box_olm_; }
    int box_olm_size() const { return box_olm_[0]->size(); }
    std::vector<std::shared_ptr<ZMatrix>>& box_olm() { return box_olm_; }
    std::shared_ptr<ZMatrix> box_olm(const int i) const { return box_olm_[i]; }

    void print_box() const;
};

}
#endif
