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

    double thresh_, extent_;
    int nbasis0_, nbasis1_, ndim_;
    int nmult_;
    double energy_ff_;

    void init();
    void insert_sp(std::vector<std::shared_ptr<const ShellPair>>);
    void insert_child(std::shared_ptr<const Box> = NULL);
    void insert_parent(std::shared_ptr<const Box> parent = NULL);
    bool is_neigh(std::shared_ptr<const Box> b, const int ws) const;
    void get_neigh(std::vector<std::shared_ptr<Box>> box, const int ws);
    void get_inter(std::vector<std::shared_ptr<Box>> box, const int ws);

    std::vector<std::complex<double>> multipole_;
    std::vector<std::complex<double>> localJ_;
    std::vector<int> offset0_, offset1_;
    void compute_M2M(std::shared_ptr<const Matrix> density);
    void sort_sp();
    std::vector<std::complex<double>> shift_multipoles(std::vector<std::complex<double>> oa, std::array<double, 3> rab) const;
    std::vector<std::complex<double>> shift_localL(std::vector<std::complex<double>> mr, std::array<double, 3> rb) const;
    std::vector<std::complex<double>> shift_localM(std::vector<std::complex<double>> olm, std::array<double, 3> r12) const;
    void compute_M2L();
    void compute_L2L();
    void compute_energy_ff();
    double compute_exact_energy_ff(std::shared_ptr<const Matrix> density) const; //debug
    std::shared_ptr<const ZMatrix> compute_node_energy(std::shared_ptr<const Matrix> density, std::vector<double> max_den, const double schwarz_thresh = 0.0) const;


  public:
    Box(int n, int id, std::array<int, 3> v, const int lmax = 10, std::vector<std::shared_ptr<const ShellPair>> sp = std::vector<std::shared_ptr<const ShellPair>>()) : rank_(n), boxid_(id), lmax_(lmax), tvec_(v), sp_(sp) { }

    ~Box() { }

    std::array<double, 3> centre() const { return centre_; }
    double centre(const int i) const { return centre_[i]; }

    int rank() const { return rank_; }
    int boxid() const { return boxid_; }
    std::array<int, 3> tvec() const { return tvec_; }

    int nsp() const { return sp_.size(); }
    int nchild() const { return child_.size(); }
    int nneigh() const { return neigh_.size(); }
    int ninter() const { return inter_.size(); }
    double extent() const { return extent_; }
    int nbasis0() const { return nbasis0_; }
    int nbasis1() const { return nbasis1_; }
    int ndim() const { return ndim_; }
    std::vector<int> offset0() const { return offset0_; }
    std::vector<int> offset1() const { return offset1_; }

    std::shared_ptr<const Box> parent() const { return parent_; }
    std::shared_ptr<const Box> child(const int i) const { return child_[i].lock(); }
    std::vector<std::shared_ptr<const Box>> neigh() const { return neigh_; }
    std::vector<std::shared_ptr<const Box>> interaction_list() const { return inter_; }
    std::vector<std::shared_ptr<const ShellPair>> sp() const { return sp_; }
    std::shared_ptr<const ShellPair> sp(const int i) const { return sp_[i]; }
    std::array<std::shared_ptr<const Shell>, 2> shells(const int i) const { return sp_[i]->shells(); }

    std::vector<std::complex<double>> multipole() const { return multipole_; }
    std::vector<std::complex<double>> localJ() const { return localJ_; }
    double energy_ff() const { return energy_ff_; }


    void print_box() const;
};

}
#endif
