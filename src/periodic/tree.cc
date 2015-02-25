//
// BAGEL - Parallel electron correlation program.
// Filename: tree.cc
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


#include <src/periodic/tree.h>
#include <src/periodic/node.h>

using namespace bagel;
using namespace std;

Tree::Tree(shared_ptr<const Molecule> mol, const int maxht) : mol_(mol), max_height_(maxht) {

  nvertex_ = mol->natom();
  position_ = mol->charge_center();
  cout << "Charge centre: " << position_[0] << "  " << position_[1] << "  " << position_[2] << endl;
  coordinates_.resize(nvertex_);
  for (int i = 0; i != nvertex_; ++i)
    coordinates_[i] = mol->atoms(i)->position();

  init();
  build_tree();
}


void Tree::init() {

  for (int i = 0; i != nvertex_; ++i) {
    coordinates_[i][0] -= position_[0];
    coordinates_[i][1] -= position_[1];
    coordinates_[i][2] -= position_[2];
  }

  get_particle_key();
//  for (int i = 0; i != nvertex_; ++i)
//    cout << particle_keys_[i] << endl;
//  cout << " Now sort " << endl;
  keysort();
  for (int i = 0; i != nvertex_; ++i)
    cout << particle_keys_[i] << endl;
}


void Tree::build_tree() {

  nnode_ = 1;
  nodes_.resize(nnode_);
  nodes_[0] = make_shared<Node>();

  bitset<nbit__>  current_key;
  for (int i = max_height_ - 1; i >= 0; --i) { /* top down */
    const int level = max_height_ - i;

    const unsigned int shift = i * 3;
    bitset<nbit__> key;

    for (int n = 0; n != nvertex_; ++n) {
      key = (leaves_[n]->key() >> shift);

      if (key != current_key) { /* insert node */
        current_key = key;
        nodes_.resize(nnode_ + 1);
        const int depth = max_height_ - i;

        if (level == 1) {
          nodes_[nnode_] = make_shared<Node>(key, depth, nodes_[0]);
          (nodes_[nnode_])->insert_vertex(leaves_[n]);
          nodes_[0]->insert_child(nodes_[nnode_]);
        } else {
          bitset<nbit__> parent_key;
          parent_key  = parent_key | (key >> 3);
          bool parent_found = false;
          for (int j = 0; j != nnode_; ++j) {
            if (parent_key == nodes_[j]->key()) {
              nodes_[nnode_] = make_shared<Node>(key, depth, nodes_[j]);
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
      }

    } // end of vertex loop
  } // end if bit loop


  print_tree_xyz();
  for (int i = 1; i != nnode_; ++i) {
    nodes_[i]->init();
    //cout << i << "  " << nodes_[i]->position(0) << "  " << nodes_[i]->position(1) << "   " << nodes_[i]->position(2) << endl;
    //cout << "Node " << i << "   has    " << nodes_[i]->nchild() << " children and is_leaf_ is " << nodes_[i]->is_leaf() << endl;
  }
}


void Tree::get_particle_key() {

  particle_keys_.resize(nvertex_);
  box_length_ = 0.0;

  const unsigned int nbitx = (nbit__ - 1) / 3;

  int iat = 0;
  for (auto& pos : coordinates_) {
    bitset<nbit__> key;
    key[nbit__-1] = 1;
    bitset<nbitx> binx(pos[0]);
    bitset<nbitx> biny(pos[1]);
    bitset<nbitx> binz(pos[2]);
    //cout << binx << " *** " << biny << " *** "  << binz << endl;
    for (int i = 0; i != 21; ++i) {
      key[i * 3 + 0] = binx[i];
      key[i * 3 + 1] = biny[i];
      key[i * 3 + 2] = binz[i];
    }
    particle_keys_[iat] = key;
    ++iat;

    const double tmp = max(pos[0], max(pos[1], pos[2]));
    if (tmp > box_length_)
      box_length_ = 2 * tmp;
  }
}


void Tree::keysort() {

  vector<int> id(nvertex_), id0(nvertex_), id1(nvertex_);
  iota(id.begin(), id.begin()+nvertex_, 0);

  for (int i = 0; i != nbit__ - 1; ++i) { // Morton order
    vector<bitset<nbit__>> key0(nvertex_), key1(nvertex_);
    int n0 = 0, n1 = 0;
    for (int n = 0; n != nvertex_; ++n) {
      if (particle_keys_[n][i] == 0) {
        id0[n0] = id[n];
        key0[n0++] = particle_keys_[n];
      } else {
        id1[n1] = id[n];
        key1[n1++] = particle_keys_[n];
      }
    }
    assert(n0 + n1 == nvertex_);
    move(key0.begin(), key0.begin()+n0, particle_keys_.begin());
    move(key1.begin(), key1.begin()+n1, particle_keys_.begin()+n0);
    move(id0.begin(), id0.begin()+n0, id.begin());
    move(id1.begin(), id1.begin()+n1, id.begin()+n0);
  }

  leaves_.resize(nvertex_);
  for (int n = 0; n != nvertex_; ++n) {
    auto leaf = make_shared<Vertex>(particle_keys_[n], mol_->atoms(id[n]));
    leaves_[n] = leaf;
  }
}


void Tree::print_tree_xyz() const { // to visualize with VMD, but not enough atoms for higher level!

  int current_level = 0;
  int node = 0;
  for (int i = 1; i != nnode_; ++i) {
    if (nodes_[i]->depth() != current_level) {
      ++current_level;
      cout << endl;
      cout << mol_->natom() << endl;
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
    for (int j = 0; j != nodes_[i]->nbody(); ++j)
      cout << setw(5) << symbol << setprecision(5) << setw(10) << nodes_[i]->bodies(j)->position(0) << "   "
                                                   << setw(10) << nodes_[i]->bodies(j)->position(1) << "   "
                                                   << setw(10) << nodes_[i]->bodies(j)->position(2) << endl;
  }
}
