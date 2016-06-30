//
// BAGEL - Parallel electron correlation program.
// Filename: tree.h
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


#ifndef __SRC_PERIODIC_TREE_H
#define __SRC_PERIODIC_TREE_H

#include <set>
#include <src/wfn/geometry.h>
#include <src/periodic/node.h>
#include <src/periodic/vertex.h>
#include <src/util/atommap.h>
#include <src/mat1e/matrix1e.h>

namespace bagel {

class Node;
class Tree {
  protected:
    std::shared_ptr<const Geometry>geom_;
    int max_height_;
    bool do_contraction_;
    int nvertex_;
    int nbasis_;
    std::vector<std::array<double, 3>> coordinates_;
    std::array<double, 3> position_;

    std::vector<std::bitset<nbit__>> particle_keys_;
    std::vector<std::shared_ptr<const Vertex>> leaves_;
    std::vector<int> ordering_, shell_id_;
    int nnode_;
    std::vector<std::shared_ptr<Node>> nodes_;
    int height_;
    // to define well-separated distributions
    double thresh_;
    int ws_;

    // vertex contraction
    std::vector<std::shared_ptr<const AtomGroup>> atomgroup_;
    void contract_vertex();


    void init();
    void build_tree();
    void get_particle_key(); // a place holder and (nbit__-1)/3 per coordinate
    void keysort();

    std::shared_ptr<const ZMatrix> compute_interactions(const int lmax, std::shared_ptr<const Matrix> density,
                             const std::vector<double> schwarz = std::vector<double>(), const double schwarz_thresh = 0.0) const;
    std::vector<std::complex<double>> get_mlm(const int lmax, std::array<double, 3> r01, std::vector<std::complex<double>> omega0) const;
    std::shared_ptr<const ZMatrix> compute_JK(std::shared_ptr<const Matrix> density) const;

  public:
    Tree(std::shared_ptr<const Geometry> geom, const int max_height = (nbit__ - 1)/3, const bool do_contract = false,
         const double thresh = PRIM_SCREEN_THRESH, const int ws = 2);
    ~Tree() { }

    void init_fmm(const int lmax, const bool dodf, const std::string auxfile) const;
    std::shared_ptr<const ZMatrix> fmm(const int lmax, std::shared_ptr<const Matrix> density = nullptr, const bool dodf = false, const double scale = 1.0, const std::vector<double> schwarz = std::vector<double>(), const double schwarz_thresh = 0.0) const;

    void print_tree_xyz() const;
    void print_leaves() const;
};


class TreeNAI : public Matrix1e {
  protected:
    void computebatch(const std::array<std::shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, std::shared_ptr<const Molecule> mol) {
      const int dimb1 = input[0]->nbasis();
      const int dimb0 = input[1]->nbasis();
      {
        NAIBatch nai(input, mol);
        nai.compute();

        add_block(2.0, offsetb1, offsetb0, dimb1, dimb0, nai.data());
      }
    }

  public:
    TreeNAI(std::shared_ptr<const Molecule> mol) : Matrix1e(mol) {
      init(mol);
      fill_upper();
    }
};

}
#endif
