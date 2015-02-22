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

Tree::Tree(shared_ptr<const Molecule> mol, const int maxht) : max_height_(maxht) {

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

  int nnode = 1;
  nodes_.resize(nnode);
  nodes_[0] = make_shared<Node>();

  bitset<nbit__>  current_key;
  for (int i = (nbit__ - 1)/3 - 1; i >= 0; --i) { /* top down */
    const int level = (nbit__ - 1)/3 - i;

    cout << "Level " << level << endl;
    if (level <= max_height_) {
      const unsigned int shift = i * 3;
      bitset<nbit__> key;

      for (int n = 0; n != nvertex_; ++n) {
        key = (leaves_[n]->key() >> shift);

        if (key != current_key) { /* insert node */
          current_key = key;
          bitset<nbit__> parent_key;
          parent_key  = parent_key | (key >> 3);
          nodes_.resize(nnode + 1);
          const int depth = (nbit__ - 1)/3 - i;
          bool parent_found = false;
          for (int j = 0; j != nnode; ++j) {
            if (parent_key == nodes_[j]->key()) {
              nodes_[nnode] = make_shared<Node>(key, depth, nodes_[j]);
              (nodes_[nnode])->insert_vertex(leaves_[n]);
              parent_found = true;
            }
          }

          if (!parent_found)
            throw runtime_error("add a node without a parent");

          ++nnode;
        } else { /* node exists */
          (nodes_[nnode-1])->insert_vertex(leaves_[n]);
        }

      } // end of vertex loop
    } // end of max height check
  } // end if bit loop

  for (int i = 1; i != nnode; ++i) {
    cout << "Level = " << nodes_[i]->depth() << " Node " << i << "  key " << nodes_[i]->key() << endl;//" parent " << (nodes_[i]->parent())->key() << endl;
    for (int j = 0; j != nodes_[i]->nbody(); ++j)
      cout << nodes_[i]->bodies(j)->position(0) << "   " << nodes_[i]->bodies(j)->position(1) << "   " << nodes_[i]->bodies(j)->position(2) << endl;
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
    auto leaf = make_shared<Vertex>(particle_keys_[n], coordinates_[id[n]]);
    leaves_[n] = leaf;
  }
}
