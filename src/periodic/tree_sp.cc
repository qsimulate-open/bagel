//
// BAGEL - Parallel electron correlation program.
// Filename: tree_sp.cc
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


#include <src/periodic/tree_sp.h>
#include <src/util/taskqueue.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/periodic/multipolebatch.h>
#include <src/periodic/localexpansion.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;
const static Legendre plm;

TreeSP::TreeSP(shared_ptr<const Geometry> geom, const int maxht, const double thresh, const int ws)
 : geom_(geom), max_height_(maxht), thresh_(thresh), ws_(ws) {

  init();
  print_leaves();
}


void TreeSP::init() {

  assert(max_height_ <= (nbit__ - 1)/3);
  centre_ = geom_->charge_center();
  nvertex_ = geom_->nshellpair();
  coordinates_.resize(nvertex_);

  nbasis_ = geom_->nbasis();
  for (int i = 0; i != nvertex_; ++i) {
    coordinates_[i][0] = geom_->shellpair(i)->centre(0) - centre_[0];
    coordinates_[i][1] = geom_->shellpair(i)->centre(1) - centre_[1];
    coordinates_[i][2] = geom_->shellpair(i)->centre(2) - centre_[2];
  }

  ordering_.resize(nvertex_);
  iota(ordering_.begin(), ordering_.begin()+nvertex_, 0);

  get_particle_key();
  keysort();

  build_tree();
  //print_tree_xyz();
  nleaf_ = 0;
  for (int i = 0; i != nnode_; ++i)
    if (nodes_[i]->is_leaf()) ++nleaf_;
}


void TreeSP::build_tree() {

  Timer treetime;
  nnode_ = 1;
  nodes_.resize(nnode_);
  nodes_[0] = make_shared<NodeSP>();

  bitset<nbit__>  current_key;
  const int max_height = max_height_;
  for (int i = 1; i <= max_height; ++i) { /* top down */
    const int depth = i;

    const unsigned int shift = nbit__ - 1 - i * 3;

    int max_nvertex = 0;
    for (int n = 0; n != nvertex_; ++n) {
      const bitset<nbit__> key = (leaves_[n]->key() >> shift);

      if (key != current_key) { /* insert node */
        current_key = key;
        nodes_.resize(nnode_ + 1);

        if (depth == 1) {
          nodes_[nnode_] = make_shared<NodeSP>(key, 1, nodes_[0], thresh_);
          (nodes_[nnode_])->insert_vertex(leaves_[n]);
          nodes_[0]->insert_child(nodes_[nnode_]);
        } else {
          const bitset<nbit__> parent_key = key >> 3;
          bool parent_found = false;
          for (int j = 0; j != nnode_; ++j) {
            if (parent_key == nodes_[j]->key()) {
              nodes_[nnode_] = make_shared<NodeSP>(key, depth, nodes_[j], thresh_);
              (nodes_[nnode_])->insert_vertex(leaves_[n]);
              parent_found = true;
              nodes_[j]->insert_child(nodes_[nnode_]);
            }
          }

          if (!parent_found)
            throw runtime_error("add a node without a parent");
        }

        ++nnode_;
      } else { /* node exists */
        (nodes_[nnode_-1])->insert_vertex(leaves_[n]);
        if (nodes_[nnode_-1]->nvertex() > max_nvertex)
          max_nvertex = nodes_[nnode_-1]->nvertex();
      }

    } // end of vertex loop
    if (max_nvertex <= 1) {
      max_height_ = i;
      break;
    }
  } // end if bit loop
  height_ = nodes_[nnode_-1]->depth();

  for (int i = 1; i != nnode_; ++i)
    nodes_[i]->init();

  for (int i = 1; i != nnode_; ++i) {
    for (int j = 1; j != nnode_; ++j) {
      if (nodes_[j]->depth() == nodes_[i]->depth()) {
        nodes_[i]->insert_neigh(nodes_[j], false, ws_);
      } else if (nodes_[j]->depth() > nodes_[i]->depth()) {
        break;
      }
    }
  }

  cout << "    * Tree construction: " << setw(15) << setprecision(2) << treetime.tick() << endl;
}


void TreeSP::init_fmm(const int lmax, const bool dodf, const string auxfile) const {
  if (dodf && auxfile.empty())
    throw runtime_error("Do FMM with DF but no df basis provided");

  Timer fmmtime;
  // Downward pass
  int u = 0;
  for (int i = nnode_ - 1; i > 0; --i) {
    if (u++ % mpi__->size() == mpi__->rank())
      nodes_[i]->compute_multipoles(lmax);
  }
  fmmtime.tick_print("Downward pass");
}


shared_ptr<const ZMatrix> TreeSP::fmm(const int lmax, shared_ptr<const Matrix> density, const bool dodf, const double scale, const vector<double> schwarz, const double schwarz_thresh) const {

  // Upward pass
  int u = 0;
  auto out = make_shared<ZMatrix>(nbasis_, nbasis_);;
  for (int i = 1; i != nnode_; ++i) {
    if (u++ % mpi__->size() == mpi__->rank()) {
//    nodes_[i]->compute_local_expansions(density, lmax);
      if (nodes_[i]->is_leaf()) {
        nodes_[i]->compute_local_expansions(density, lmax); //////// TMP
        shared_ptr<const ZMatrix> tmp = nodes_[i]->compute_Coulomb(density, dodf, schwarz, schwarz_thresh);
        *out += *tmp;
      }
    }
  }

  // return the Coulomb matrix
  out->allreduce();

  return out;
}


void TreeSP::get_particle_key() {

  particle_keys_.resize(nvertex_);

  const unsigned int nbitx = (nbit__ - 1) / 3;

  int isp = 0;
  for (auto& pos : coordinates_) {
    bitset<nbit__> key;
    key[nbit__-1] = 1;
    bitset<nbitx> binx(pos[0]);
    bitset<nbitx> biny(pos[1]);
    bitset<nbitx> binz(pos[2]);
    for (int i = 0; i != 21; ++i) {
      key[i * 3 + 0] = binx[i];
      key[i * 3 + 1] = biny[i];
      key[i * 3 + 2] = binz[i];
    }
    particle_keys_[isp] = key;
    ++isp;
  }
}


void TreeSP::keysort() {

  vector<int> id0(nvertex_), id1(nvertex_);

  for (int i = 0; i != nbit__ - 1; ++i) { // Morton order
    vector<bitset<nbit__>> key0(nvertex_), key1(nvertex_);
    int n0 = 0, n1 = 0;
    for (int n = 0; n != nvertex_; ++n) {
      if (particle_keys_[n][i] == 0) {
        id0[n0] = ordering_[n];
        key0[n0++] = particle_keys_[n];
      } else {
        id1[n1] = ordering_[n];
        key1[n1++] = particle_keys_[n];
      }
    }
    assert(n0 + n1 == nvertex_);
    move(key0.begin(), key0.begin()+n0, particle_keys_.begin());
    move(key1.begin(), key1.begin()+n1, particle_keys_.begin()+n0);
    move(id0.begin(), id0.begin()+n0, ordering_.begin());
    move(id1.begin(), id1.begin()+n1, ordering_.begin()+n0);
  }

  leaves_.resize(nvertex_);
  for (int n = 0; n != nvertex_; ++n) {
    const int pos = ordering_[n];
    auto leaf = make_shared<VertexSP>(particle_keys_[n], geom_->shellpair(pos));
    leaves_[n] = leaf;
  }
}


void TreeSP::print_tree_xyz() const { // to visualize with VMD, but not enough atoms for higher level!

  int current_level = 0;
  int node = 0;
  for (int i = 1; i != nnode_; ++i) {
    if (nodes_[i]->depth() != current_level) {
      ++current_level;
      cout << endl;
      cout << nvertex_ << endl;
      cout << "Level " << current_level << endl;
      node = 0;
    }
    ++node;
    string symbol("");
    switch(node) {
      case 1: symbol = "H "; break;
      case 2: symbol = "He"; break;
      case 3: symbol = "Li"; break;
      case 4: symbol = "Be"; break;
      case 5: symbol = "B "; break;
      case 6: symbol = "C "; break;
      case 7: symbol = "N "; break;
      case 8: symbol = "O "; break;
      case 9: symbol = "F "; break;
      case 10: symbol = "Ne"; break;
      case 11: symbol = "Na"; break;
      case 12: symbol = "Mg"; break;
      case 13: symbol = "Al"; break;
      case 14: symbol = "Si"; break;
      case 15: symbol = "P "; break;
      case 16: symbol = "S "; break;
      case 17: symbol = "Cl"; break;
      case 18: symbol = "Ar"; break;
      case 19: symbol = "K "; break;
      case 20: symbol = "Ca"; break;
      case 21: symbol = "Sc"; break;
      case 22: symbol = "Ti"; break;
      case 23: symbol = "V "; break;
      case 24: symbol = "Cr"; break;
      case 25: symbol = "Mn"; break;
      case 26: symbol = "Fe"; break;
      case 27: symbol = "Co"; break;
      case 28: symbol = "Ni"; break;
      case 29: symbol = "Cu"; break;
      case 30: symbol = "Zn"; break;
      case 31: symbol = "Ga"; break;
      case 32: symbol = "Ge"; break;
      case 33: symbol = "As"; break;
      case 34: symbol = "Se"; break;
      case 35: symbol = "Br"; break;
      case 36: symbol = "Kr"; break;
      case 37: symbol = "Rb"; break;
      case 38: symbol = "Sr"; break;
      case 39: symbol = "Y "; break;
      case 40: symbol = "Zr"; break;
      case 41: symbol = "Nb"; break;
      case 42: symbol = "Mo"; break;
      case 43: symbol = "Tc"; break;
      case 44: symbol = "Ru"; break;
      case 45: symbol = "Rh"; break;
      case 46: symbol = "Pd"; break;
      case 47: symbol = "Ag"; break;
      case 48: symbol = "Cd"; break;
      case 49: symbol = "In"; break;
      case 50: symbol = "Sn"; break;
      case 51: symbol = "Sb"; break;
      case 52: symbol = "Te"; break;
      case 53: symbol = "I "; break;
      case 54: symbol = "Xe"; break;
      case 55: symbol = "Cs"; break;
      case 56: symbol = "Ba"; break;
      case 57: symbol = "La"; break;
      case 58: symbol = "Ce"; break;
      case 59: symbol = "Pr"; break;
      case 60: symbol = "Nd"; break;
      case 61: symbol = "Pm"; break;
      case 62: symbol = "Sm"; break;
      case 63: symbol = "Eu"; break;
      case 64: symbol = "Gd"; break;
    }
    for (int j = 0; j != nodes_[i]->nvertex(); ++j) {
      cout << setw(5) << symbol << setprecision(5) << setw(10) << nodes_[i]->vertex(j)->centre(0) << "   "
                                                   << setw(10) << nodes_[i]->vertex(j)->centre(1) << "   "
                                                   << setw(10) << nodes_[i]->vertex(j)->centre(2) << endl;
      if (nodes_[i]->depth() == max_height_) {
        cout << "*** Neighbours: " << endl;
        for (auto& neigh : nodes_[i]->neigh())
          cout << setprecision(5) << setw(10) << neigh->centre(0) << "  "
                                  << setw(10) << neigh->centre(1) << "  "
                                  << setw(10) << neigh->centre(2) << endl;
      }
    }
  }
}


void TreeSP::print_leaves() const {

  const std::string space = "       ";
  cout << "   i      Rank   Nvertex   Ninter    Nneigh   Extent " << endl;
  for (int i = 0; i != nnode_; ++i)
    if (nodes_[i]->is_leaf())
      cout << "  " << i << space << nodes_[i]->rank() << space << nodes_[i]->nvertex() << space
           << nodes_[i]->interaction_list().size() << space << nodes_[i]->nneigh() << space << nodes_[i]->extent() << endl;
  cout << "#NODES = " << nnode_ << " #LEAVES = " << nleaf_ << endl;
}
