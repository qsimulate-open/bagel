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
#include <src/util/taskqueue.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/periodic/multipolebatch.h>
#include <src/periodic/localexpansion.h>

using namespace bagel;
using namespace std;

static const AtomMap atommap_;
static const double pisq__ = pi__ * pi__;
const static Legendre plm;

Tree::Tree(shared_ptr<const Geometry> geom, const int maxht, const bool do_contract, const int lmax, const double thresh, const int ws)
 : geom_(geom), schwarz_(geom_->schwarz()), max_height_(maxht), do_contraction_(do_contract), lmax_(lmax), thresh_(thresh), ws_(ws) {

  init();
  print_leaves();
}


void Tree::init() {

  assert(max_height_ <= (nbit__ - 1)/3);
  position_ = geom_->charge_center();
  atomgroup_.resize(geom_->natom());
  coordinates_.resize(geom_->natom());

  nbasis_ = geom_->nbasis();
  shell_id_.resize(geom_->natom());
  shell_id_[0] = 0;
  for (int i = 1; i != geom_->natom(); ++i)
    shell_id_[i] = shell_id_[i-1] + geom_->atoms(i-1)->nshell();

  if (do_contraction_) {
    contract_vertex();
  } else {
    nvertex_ = geom_->natom();
    for (int i = 0; i != nvertex_; ++i) {
      atomgroup_[i] = make_shared<const AtomGroup>
                      (vector<shared_ptr<const Atom>>(1, geom_->atoms(i)), vector<int>(1, i), vector<int>(1, shell_id_[i]));
      coordinates_[i][0] = geom_->atoms(i)->position(0) - position_[0];
      coordinates_[i][1] = geom_->atoms(i)->position(1) - position_[1];
      coordinates_[i][2] = geom_->atoms(i)->position(2) - position_[2];
    }
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


void Tree::contract_vertex() {

  constexpr double pm2au = 1e-10 / au2meter__;
  constexpr double tol = 0.0;
  const set<string> svalence = {"h", "li", "f", "na", "cl", "k", "br", "i"};

  int nvertex = 0;
  vector<int> cont_vertex(geom_->natom(), -1);
  for (int i = 0; i != geom_->natom(); ++i) {
    if (cont_vertex[i] < 0) {
      const bool is_singly_valent_atom = svalence.find(geom_->atoms(i)->name()) != svalence.end();
      if (is_singly_valent_atom) {
        const double rad_i = pm2au * atommap_.cov_radius(geom_->atoms(i)->name());

        vector<shared_ptr<const Atom>> tmpatom;
        vector<int> tmporder;
        vector<int> tmpshell;
        tmpatom.insert(tmpatom.end(), geom_->atoms(i));
        tmporder.insert(tmporder.end(), i);
        tmpshell.insert(tmpshell.end(), shell_id_[i]);
        bool is_found = false;
        for (int j = 0; j != geom_->natom(); ++j) {
          if (j != i) {
            array<double, 3> v12;
            v12[0] = geom_->atoms(i)->position(0) - geom_->atoms(j)->position(0);
            v12[1] = geom_->atoms(i)->position(1) - geom_->atoms(j)->position(1);
            v12[2] = geom_->atoms(i)->position(2) - geom_->atoms(j)->position(2);
            const double bondlength = sqrt(v12[0]*v12[0] + v12[1]*v12[1] + v12[2]*v12[2]);
            const double rad_j = pm2au * atommap_.cov_radius(geom_->atoms(j)->name());
            if (bondlength <= (1.0 + tol) * (rad_i + rad_j)) {
              const int ivertex = cont_vertex[j];
              is_found = true;
              if (ivertex < 0) {
                tmpatom.insert(tmpatom.end(), geom_->atoms(j));
                tmporder.insert(tmporder.end(), j);
                tmpshell.insert(tmpshell.end(), shell_id_[j]);
                cont_vertex[i] = nvertex;
                cont_vertex[j] = nvertex;

                atomgroup_[nvertex] = make_shared<const AtomGroup>(tmpatom, tmporder, tmpshell);
                ++nvertex;
              } else {
                cont_vertex[i] = ivertex;
                const vector<shared_ptr<const Atom>> atom = atomgroup_[ivertex]->atoms();
                const vector<int> order = atomgroup_[ivertex]->order_in_geom();
                const vector<int> shell = atomgroup_[ivertex]->ishell();
                tmpatom.insert(tmpatom.end(), atom.begin(), atom.end());
                tmporder.insert(tmporder.end(), order.begin(), order.end());
                tmpshell.insert(tmpshell.end(), shell.begin(), shell.end());

                atomgroup_[ivertex] = make_shared<const AtomGroup>(tmpatom, tmporder, tmpshell);
              }
              break;
            }
          }
        }
        if (!is_found)
          throw runtime_error("Found a radical");

      } else {
        atomgroup_[nvertex] = make_shared<const AtomGroup>
                              (vector<shared_ptr<const Atom>>(1, geom_->atoms(i)), vector<int>(1, i), vector<int>(1, shell_id_[i]));
        cont_vertex[i] = nvertex;
        ++nvertex;
      }
    }
  }

  nvertex_ = nvertex;
  atomgroup_.resize(nvertex_);
  coordinates_.resize(nvertex_);
  for (int i = 0; i != nvertex_; ++i) {
    coordinates_[i][0] = atomgroup_[i]->position(0) - position_[0];
    coordinates_[i][1] = atomgroup_[i]->position(1) - position_[1];
    coordinates_[i][2] = atomgroup_[i]->position(2) - position_[2];
//    for (int iat = 0; iat != atomgroup_[i]->atoms().size(); ++iat)
//      cout << atomgroup_[i]->atoms(iat)->name() << atomgroup_[i]->order_in_geom(iat) << " ";
//    cout << endl;
  }
}


void Tree::build_tree() {

  Timer treetime;
  nnode_ = 1;
  nodes_.resize(nnode_);
  nodes_[0] = make_shared<Node>();

  bitset<nbit__>  current_key;
  const int max_height = max_height_;
  for (int i = 1; i <= max_height; ++i) { /* top down */
    const int depth = i;

    const unsigned int shift = nbit__ - 1 - i * 3;

    int max_nbody = 0;
    for (int n = 0; n != nvertex_; ++n) {
      const bitset<nbit__> key = (leaves_[n]->key() >> shift);

      if (key != current_key) { /* insert node */
        current_key = key;
        nodes_.resize(nnode_ + 1);

        if (depth == 1) {
          nodes_[nnode_] = make_shared<Node>(key, 1, nodes_[0], thresh_);
          (nodes_[nnode_])->insert_vertex(leaves_[n]);
          nodes_[0]->insert_child(nodes_[nnode_]);
        } else {
          const bitset<nbit__> parent_key = key >> 3;
          bool parent_found = false;
          for (int j = 0; j != nnode_; ++j) {
            if (parent_key == nodes_[j]->key()) {
              nodes_[nnode_] = make_shared<Node>(key, depth, nodes_[j], thresh_);
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
        if (nodes_[nnode_-1]->nbody() > max_nbody)
          max_nbody = nodes_[nnode_-1]->nbody();
      }

    } // end of vertex loop
    if (max_nbody <= 1) {
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
        nodes_[i]->insert_neighbour(nodes_[j], false, ws_);
      } else if (nodes_[j]->depth() > nodes_[i]->depth()) {
        break;
      }
    }
  }

  cout << "    * Tree construction: " << setw(15) << setprecision(2) << treetime.tick() << endl;
}


void Tree::init_fmm(const bool dodf, const string auxfile) const {

  if (lmax_ < 0) return;
  if (dodf && auxfile.empty())
    throw runtime_error("Do FMM with DF but no df basis provided");

  Timer fmmtime;
#if 1
  // Downward pass
  int u = 0;
  for (int i = nnode_ - 1; i > 0; --i) {
    if (u++ % mpi__->size() == mpi__->rank()) {
      nodes_[i]->compute_multipoles(lmax_);
      if (dodf && nodes_[i]->is_leaf())
        nodes_[i]->form_df(auxfile);
    }
  }
#endif
  fmmtime.tick_print("Downward pass");
}


shared_ptr<const ZMatrix> Tree::fmm(shared_ptr<const Matrix> density, const bool dodf, const double scale, const double schwarz_thresh) const {

#if 0
  ///////// DEBUG ///////////
  shared_ptr<const ZMatrix> out = compute_interactions(lmax_, density, schwarz_, schwarz_thresh);
  /////////// END OF DEBUG ///////////
#endif

  #if 1
  // Upward pass
  vector<int> offsets;
  for (int n = 0; n != geom_->natom(); ++n) {
    const vector<int> tmpoff = geom_->offset(n);
    offsets.insert(offsets.end(), tmpoff.begin(), tmpoff.end());
  }

  int u = 0;
  auto out = make_shared<ZMatrix>(nbasis_, nbasis_);;
  int nint = 0;
  for (int i = 1; i != nnode_; ++i) {
    if (u++ % mpi__->size() == mpi__->rank()) {
//    nodes_[i]->compute_local_expansions(density, lmax_, offsets, scale);
      if (nodes_[i]->is_leaf()) {
        nodes_[i]->compute_local_expansions(density, lmax_, offsets, scale); //////// TMP
        shared_ptr<const ZMatrix> tmp = nodes_[i]->compute_Coulomb(nbasis_, density, offsets, dodf, scale, schwarz_, schwarz_thresh);
        *out += *tmp;
        nint += nodes_[i]->shellquads_.size();
      }
    }
  }

  // return the Coulomb matrix
  out->allreduce();
//  shared_ptr<const ZMatrix> nf = compute_JK(density, nint);
//  *out += *nf;

  #endif
  return out;
}


void Tree::get_particle_key() {

  particle_keys_.resize(nvertex_);

  const unsigned int nbitx = (nbit__ - 1) / 3;

  int iat = 0;
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
    particle_keys_[iat] = key;
    ++iat;
  }
}


void Tree::keysort() {

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
    auto leaf = make_shared<Vertex>(particle_keys_[n], atomgroup_[pos]);
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
    for (int j = 0; j != nodes_[i]->nbody(); ++j) {
      cout << setw(5) << symbol << setprecision(5) << setw(10) << nodes_[i]->bodies(j)->position(0) << "   "
                                                   << setw(10) << nodes_[i]->bodies(j)->position(1) << "   "
                                                   << setw(10) << nodes_[i]->bodies(j)->position(2) << endl;
      if (nodes_[i]->depth() == max_height_) {
        cout << "*** Neighbours: " << endl;
        for (auto& neigh : nodes_[i]->neighbour())
          cout << setprecision(5) << setw(10) << neigh->position(0) << "  "
                                  << setw(10) << neigh->position(1) << "  "
                                  << setw(10) << neigh->position(2) << endl;
      }
    }
  }
}


void Tree::print_leaves() const {

  const std::string space = "       ";
  cout << "   i      Rank   Nbody   Ninter    Nneigh   Extent " << endl;
  for (int i = 0; i != nnode_; ++i)
    if (nodes_[i]->is_leaf())
      cout << "  " << i << space << nodes_[i]->rank() << space << nodes_[i]->nbody() << space
           << nodes_[i]->interaction_list().size() << space << nodes_[i]->nneighbour() << space << nodes_[i]->extent() << endl;
  cout << "#NODES = " << nnode_ << " #LEAVES = " << nleaf_ << endl;
}


shared_ptr<const ZMatrix> Tree::compute_interactions(shared_ptr<const Matrix> density, const double schwarz_thresh) const {

  const int nspairs = geom_->nshellpair();
  const int nsh = sqrt(nspairs);
  assert(nsh*nsh == nspairs);
  auto out = make_shared<ZMatrix>(nbasis_, nbasis_);
  if (!density) return out;

  const double* density_data = density->data();

  vector<double> max_density(nspairs);
  for (int i01 = 0; i01 != nspairs; ++i01) {
    shared_ptr<const Shell> sh0 = geom_->shellpair(i01)->shell(0);
    const int offset0 = geom_->shellpair(i01)->offset(0);
    const int i0 = geom_->shellpair(i01)->shell_ind(0);
    const int size0 = sh0->nbasis();

    shared_ptr<const Shell> sh1 = geom_->shellpair(i01)->shell(1);
    const int offset1 = geom_->shellpair(i01)->offset(1);
    const int i1 = geom_->shellpair(i01)->shell_ind(1);
    const int size1 = sh0->nbasis();

    double denmax = 0.0;
    for (int i0 = offset0; i0 != offset0 + size0; ++i0) {
      const int i0n = i0 * density->ndim();
      for (int i1 = offset1; i1 != offset1 + size1; ++i1)
        denmax = max(denmax, fabs(density_data[i0n + i1]));
    }
    max_density[i01] = denmax;
  }

 #if 1
  const int shift = sizeof(int) * 4;
  const int nmult = (lmax_+1)*(lmax_+1);
  shared_ptr<Petite> plist = geom_->plist();;
  int nzero = 0;
  int nint = 0;
  for (int i01 = 0; i01 != nspairs; ++i01) {
//    if (!plist->in_p2(i01)) continue;
    if (schwarz_[i01] < 1e-15) continue;
    const double density_01 = max_density[i01] * 4.0;

    shared_ptr<const Shell> sh0 = geom_->shellpair(i01)->shell(0);
    const int i0 = geom_->shellpair(i01)->shell_ind(0);
    const int offset0 = geom_->shellpair(i01)->offset(0);
    const int size0 = sh0->nbasis();
//    if (!plist->in_p1(i0)) continue;

    shared_ptr<const Shell> sh1 = geom_->shellpair(i01)->shell(1);
    const int offset1 = geom_->shellpair(i01)->offset(1);
    const int i1 = geom_->shellpair(i01)->shell_ind(1);
    const int size1 = sh1->nbasis();

    assert(i01 == i0 * nsh + i1);
    if (i1 < i0) continue;

    for (int i23 = i01; i23 != nspairs; ++i23) {
      if (schwarz_[i23] < 1e-15) continue;
      const double density_23 = max_density[i23] * 4.0;

      shared_ptr<const Shell> sh2 = geom_->shellpair(i23)->shell(0);
      const int i2 = geom_->shellpair(i23)->shell_ind(0);
      const int offset2 = geom_->shellpair(i23)->offset(0);
      const int size2 = sh2->nbasis();
      if (i2 < i0) continue;

      shared_ptr<const Shell> sh3 = geom_->shellpair(i23)->shell(1);
      const int offset3 = geom_->shellpair(i23)->offset(1);
      const int i3 = geom_->shellpair(i23)->shell_ind(1);
      const int size3 = sh3->nbasis();

      assert(i23 == i2 * nsh + i3);
      if (i3 < i2) continue;
//      int ijkl = plist->in_p4(i01, i23, i0, i1, i2, i3);
//      if (ijkl == 0) continue;

      const double density_02 = max_density[i0 * nsh + i2];
      const double density_03 = max_density[i0 * nsh + i3];
      const double density_12 = max_density[i1 * nsh + i2];
      const double density_13 = max_density[i1 * nsh + i3];

      const double mulfactor = max(max(max(density_01, density_02),
                                       max(density_03, density_12)),
                                       max(density_13, density_23));
      const double integral_bound = mulfactor * schwarz_[i01] * schwarz_[i23];
      const bool skip_schwarz = integral_bound < schwarz_thresh;
      if (skip_schwarz) continue;

      array<shared_ptr<const Shell>,4> input = {{sh3, sh2, sh1, sh0}};
      const bool is_neigh = geom_->shellpair(i01)->is_neighbour(geom_->shellpair(i23), ws_);

      if (!is_neigh) {
        // multipoles
        #if 0
        vector<complex<double>> omega01(nmult), omega23(nmult);
        vector<shared_ptr<const ZMatrix>> qlm01 = geom_->shellpair(i01)->multipoles();
        vector<shared_ptr<const ZMatrix>> qlm23 = geom_->shellpair(i23)->multipoles();

        auto den01 = density->get_submatrix(offset1, offset0, size1, size0);
        auto den23 = density->get_submatrix(offset3, offset2, size3, size2);
        for (int i = 0; i != nmult; ++i) {
          omega01[i] = qlm01[i]->get_real_part()->dot_product(*den01);
          omega23[i] = qlm23[i]->get_real_part()->dot_product(*den23);
        }
        array<double, 3> rvec0123, rvec2301;
        for (int i = 0; i != 3; ++i) {
          rvec0123[i] = geom_->shellpair(i23)->centre(i) - geom_->shellpair(i01)->centre(i);
          rvec2301[i] = -rvec0123[i];
        }
        vector<complex<double>> mlm01 = get_mlm(rvec0123, omega01);
        vector<complex<double>> mlm23 = get_mlm(rvec2301, omega23);
        for (int j2 = 0; j2 != size2; ++j2)
          for (int j3 = 0; j3 != size3; ++j3)
              for (int i = 0; i != nmult; ++i)
                out->element(j3+offset3, j2+offset2) += qlm23[i]->element(j3, j2) * mlm01[i];
        for (int j0 = 0; j0 != size0; ++j0)
          for (int j1 = 0; j1 != size1; ++j1)
              for (int i = 0; i != nmult; ++i)
                out->element(j1+offset1, j0+offset0) += qlm01[i]->element(j1, j0) * mlm23[i];
        #endif
        continue;
      } else {
        // integrate
        ++nint;
        ERIBatch eribatch(input, mulfactor);
        eribatch.compute();
        const double* eridata = eribatch.data();

        for (int j0 = offset0; j0 != offset0 + size0; ++j0) {
          const int j0n = j0 * density->ndim();
          for (int j1 = offset1; j1 != offset1 + size1; ++j1) {
            if (j1 < j0) {
              eridata += size2 * size3;
              continue;
            }

            const int j1n = j1 * density->ndim();
            const double scale_01 = (j0 == j1) ? 0.5 : 1.0;
            const unsigned int nj01 = (j0 << shift) + j1;

            for (int j2 = offset2; j2 != offset2 + size2; ++j2) {
              const int j2n = j2 * density->ndim();

              const int maxj0j2 = max(j0, j2);
              const int minj0j2 = min(j0, j2);

              const int maxj1j2 = max(j1, j2);
              const int minj1j2 = min(j1, j2);

              for (int j3 = offset3; j3 != offset3 + size3; ++j3, ++eridata) {
                if (j3 < j2) continue;
                const unsigned int nj23 = (j2 << shift) + j3;
                if ((nj23 < nj01) && (i01 == i23)) continue;

                const int j3n = j3 * density->ndim();
                const double scale_23 = (j2 == j3) ? 0.5 : 1.0;
                const double scale = ((nj01 == nj23) ? 0.25 : 0.5);
                //const double scale = ((nj01 == nj23) ? 0.25 : 0.5) * static_cast<double>(ijkl);

                const int maxj1j3 = max(j1, j3);
                const int minj1j3 = min(j1, j3);
                const double eri = *eridata;

                const double intval = eri * scale_01 * scale_23 * scale;
                out->element(j1, j0) += density_data[j2n + j3] * intval * 4.0;
                out->element(j3, j2) += density_data[j0n + j1] * intval * 4.0;
                out->element(j3, j0) -= density_data[j2n + j1] * intval;

                out->element(maxj1j2, minj1j2) -= density_data[j0n + j3] * intval;
                out->element(maxj0j2, minj0j2) -= density_data[j1n + j3] * intval;
                out->element(maxj1j3, minj1j3) -= density_data[j0n + j2] * intval;
              }
            }
          }
        }
      }

    }
  }
  cout << "nint in compute_interactions = " << nint << endl;
  for (int i = 0; i != density->ndim(); ++i) out->element(i, i) *= 2.0;
  out->fill_upper();
  #endif

  return out;
}


vector<complex<double>> Tree::get_mlm(array<double, 3> r01, vector<complex<double>> omega0) const {

  const int nmult = (lmax_ + 1) * (lmax_ + 1);
  assert(omega0.size() == nmult);
  const double r = sqrt(r01[0]*r01[0] + r01[1]*r01[1] + r01[2]*r01[2]);
  const double ctheta = (r > numerical_zero__) ? r01[2]/r : 0.0;
  const double phi = atan2(r01[1], r01[0]);

  vector<complex<double>> out(nmult);

  int i1 = 0;
  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m, ++i1) {


      out[i1] = 0.0;
      int i2 = 0;
      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k, ++i2) {

          const int a = l + j;
          const int b = m - l + k - j;

          double prefactor = plm.compute(a, abs(b), ctheta) / pow(r, a + 1);
          double ft = 1.0;
          for (int i = 1; i <= a - abs(b); ++i) {
            prefactor *= ft;
            ++ft;
          }
          const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (-1.0 * prefactor * cos(abs(b) * phi));
          const double imag = prefactor * sin(abs(b) * phi);
          const complex<double> coeff(real, imag);

          out[i1] += omega0[i2] * coeff;
        } //jk
      }
    } //lm
  }

  return out;
}


shared_ptr<const ZMatrix> Tree::compute_JK(shared_ptr<const Matrix> density, const int nint) const {

  auto out = make_shared<ZMatrix>(nbasis_, nbasis_);
  if (!density) return out;
  vector<shared_ptr<const Shell>> shells;
  for (auto& atom : geom_->atoms()) {
    const vector<shared_ptr<const Shell>> tmpsh = atom->shells();
    shells.insert(shells.end(), tmpsh.begin(), tmpsh.end());
  }

  TaskQueue<function<void(void)>> tasks(nleaf_);
/////  TaskQueue<function<void(void)>> tasks(nint);
  for (int i = 0; i != nnode_; ++i) {
    if (!nodes_[i]->is_leaf()) continue;
    tasks.emplace_back(
      [this, i, &out, &density, shells] () {
    int ish = 0;
    for (auto& sh : nodes_[i]->shellquads_) {
      shared_ptr<Petite> plist = geom_->plist();;
      const int nsh = shells.size();
      if (!plist->in_p1(sh[0])) continue;
      const int i01 = sh[0] * nsh + sh[1];
      const int i23 = sh[2] * nsh + sh[3];
      if (!plist->in_p2(i01)) continue;
      int ijkl = plist->in_p4(i01, i23, sh[0], sh[1], sh[2], sh[3]);
      if (ijkl == 0) continue;

      shared_ptr<const Shell> b0 = shells[sh[0]];
      shared_ptr<const Shell> b1 = shells[sh[1]];
      shared_ptr<const Shell> b2 = shells[sh[2]];
      shared_ptr<const Shell> b3 = shells[sh[3]];
      array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
      array<size_t, 4> offsets = nodes_[i]->int_offsets_[ish];
      ++ish;

      const size_t b0offset = offsets[0];
      const size_t b1offset = offsets[1];
      const size_t b2offset = offsets[2];
      const size_t b3offset = offsets[3];
      if ((b0offset + input[3]->nbasis() + b1offset + input[2]->nbasis()) < (b2offset + b3offset)) continue;
      if ((b0offset + input[3]->nbasis() + b2offset + input[1]->nbasis()) < (b1offset + b3offset)) continue;

/////      tasks.emplace_back(
/////        [this, &out, &density, input, offsets, ijkl] () {
          const double* density_data = density->data();
          ERIBatch eribatch(input, 0.0);
          eribatch.compute();
          const double* eridata = eribatch.data();


          for (int j0 = b0offset; j0 != b0offset + input[3]->nbasis(); ++j0) {
            const int j0n = j0 * density->ndim();
            for (int j1 = b1offset; j1 != b1offset + input[2]->nbasis(); ++j1) {
              const int j1n = j1 * density->ndim();
              const double scale_01 = (j0 == j1) ? 1.0 : 2.0;
              for (int j2 = b2offset; j2 != b2offset + input[1]->nbasis(); ++j2) {
                const int j2n = j2 * density->ndim();
                for (int j3 = b3offset; j3 != b3offset + input[0]->nbasis(); ++j3, ++eridata) {
                  const int j3n = j3 * density->ndim();
                  double eri = *eridata * static_cast<double>(ijkl) * 0.5;

                  if (j0 + j1 <  j2 + j3) continue;
                  if (j0 + j2 <  j1 + j3) continue;
                  if (j0 + j1 == j2 + j3) eri *= 0.5;
                  if (j0 + j2 == j1 + j3) eri *= 0.5;
                  const double eri2 = eri * 2.0;
                  out->element(j1, j0) += density_data[j3n + j2] * eri2;
                  out->element(j3, j0) -= density_data[j1n + j2] * eri;

                  out->element(j3, j2) += density_data[j1n + j0] * eri2;
                  out->element(j1, j2) -= density_data[j3n + j0] * eri;

                  out->element(j0, j1) += density_data[j2n + j3] * eri2;
                  out->element(j0, j3) -= density_data[j2n + j1] * eri;

                  out->element(j2, j3) += density_data[j0n + j1] * eri2;
                  out->element(j2, j1) -= density_data[j0n + j3] * eri;
                }
              }
            }
          }
//////        }
//////      );
    }
        }
      ); // end task
  }
  tasks.compute();
  out->allreduce();

  return out;

}
