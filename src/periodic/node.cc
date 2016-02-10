//
// BAGEL - Parallel electron correlation program.
// Filename: node.cc
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


#include <src/periodic/node.h>
#include <src/periodic/multipolebatch.h>
#include <src/integral/os/overlapbatch.h>
#include <src/periodic/localexpansion.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;

Node::Node(const bitset<nbit__> key, const int depth, shared_ptr<const Node> parent, const double thresh)
 : key_(key), depth_(depth), parent_(parent), thresh_(thresh) {
  if (depth == 0)
    key_[0] = 1;

  is_leaf_ = false;
  is_complete_ = false;
  nchild_  = 0;
  nbody_   = 0;
  nneighbour_ = 0;
  ninter_ = 0;
  is_same_as_parent_ = false;
  rank_ = 0;
}


void Node::insert_vertex(shared_ptr<const Vertex> p) {

  assert(!is_complete_);
  bodies_.resize(nbody_ + 1);
  bodies_[nbody_] = p;

  ++nbody_;
}


void Node::insert_child(shared_ptr<const Node> child) {

  if (child) {
    assert(!is_leaf_);
    children_.resize(nchild_ + 1);
    children_[nchild_] = child;
    ++nchild_;
  }
}


void Node::init() {

  is_complete_ = true;
  iself_ = -1;
  if (nchild_ == 0) is_leaf_ = true;

  if (!is_leaf_)
    child_local_expansion_.resize(nchild_);

  position_ = {{0.0, 0.0, 0.0}};
  nbasis_ = 0;
  double sum = 0.0;
  for (auto& body : bodies_) {
    for (auto& atom : body->atoms()) {
      position_[0] += atom->atom_charge() * atom->position(0);
      position_[1] += atom->atom_charge() * atom->position(1);
      position_[2] += atom->atom_charge() * atom->position(2);
      sum += atom->atom_charge();
      nbasis_ += atom->nbasis();
    }
  }
  position_[0] /= sum;
  position_[1] /= sum;
  position_[2] /= sum;

  array<double, 3> r12;
  r12[0] = position_[0] - parent_->position(0);
  r12[1] = position_[1] - parent_->position(1);
  r12[2] = position_[2] - parent_->position(2);

  const double r = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
  if (r <= numerical_zero__) {
    is_same_as_parent_ = true;
    rank_ = parent_->rank();
  } else {
    rank_ = parent_->rank() + 1;
  }

  compute_extent(thresh_);
}


void Node::compute_extent(const double thresh) {

  vector<shared_ptr<const Shell>> shells;
  for (auto& body : bodies_)
    for (auto& atom : body->atoms())
      shells.insert(shells.end(), atom->shells().begin(), atom->shells().end());

  const array<double, 3> centre = position_; // use charge centre for now

  extent_ = 0.0;
  for (auto& ish : shells)
    for (auto& jsh : shells) {
      const vector<double> exp0 = ish->exponents();
      const vector<double> exp1 = jsh->exponents();
      array<double, 3> AB;
      AB[0] = ish->position(0) - jsh->position(0);
      AB[1] = ish->position(1) - jsh->position(1);
      AB[2] = ish->position(2) - jsh->position(2);
      const double rsq = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
      const double lnthresh = log(thresh);

      for (auto& expi0 : exp0) {
        for (auto& expi1 : exp1) {
          const double cxp_inv = 1.0 / (expi0 + expi1);
          const double expi01 = expi0 * expi1;
          const double lda_kl = sqrt((- lnthresh - expi01 * rsq * cxp_inv + 0.75 * log(4.0 * expi01 / pisq__)) * cxp_inv);

          array<double, 3> tmp;
          tmp[0] = (ish->position(0) * expi0 + jsh->position(0) * expi1) * cxp_inv - centre[0];
          tmp[1] = (ish->position(1) * expi0 + jsh->position(1) * expi1) * cxp_inv - centre[1];
          tmp[2] = (ish->position(2) * expi0 + jsh->position(2) * expi1) * cxp_inv - centre[2];

          const double extent0 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]) + lda_kl;
          if (extent0 > extent_) extent_ = extent0;
        }
      }
    }
}


void Node::insert_neighbour(shared_ptr<const Node> neigh, const bool is_neighbour, const int ws) {

  assert(neigh->depth() == depth_);
  if (!is_neighbour) {
      array<double, 3> r12;
      r12[0] = position_[0] - neigh->position(0);
      r12[1] = position_[1] - neigh->position(1);
      r12[2] = position_[2] - neigh->position(2);
      const double r = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      if (r <= (1.0 + ws) * (extent_ + neigh->extent())) {
        neighbour_.resize(nneighbour_ + 1);
        neighbour_[nneighbour_] = neigh;
        if (r < numerical_zero__)
          iself_ = nneighbour_;
        ++nneighbour_;
      } else {
        interaction_list_.resize(ninter_ + 1);
        interaction_list_[ninter_] = neigh;
        ++ninter_;
      }
#if 0
    } else {
      interaction_list_.resize(ninter_ + 1);
      interaction_list_[ninter_] = neigh;
      ++ninter_;
    }
#endif
  } else {
    neighbour_.resize(nneighbour_ + 1);
    neighbour_[nneighbour_] = neigh;
    ++nneighbour_;
  }
}


void Node::make_interaction_list(const int ws) {

  ninter_ = 0;
  for (auto& pneighbour : parent_->neighbour()) {
    for (int n = 0; n != pneighbour->nchild(); ++n) {
      shared_ptr<const Node> child = pneighbour->children(n);
      assert(child->depth() == depth_);
      bool is_close = false;
      for (int i = 0; i != nneighbour_; ++i) {
        if (child->key() == neighbour_[i]->key()) {
          is_close = true;
          break;
        }
      }
      if (!is_close) {
        interaction_list_.resize(ninter_ + 1);
        interaction_list_[ninter_] = child;
        ++ninter_;
      }
    }
  }
}


array<double, 3> Node::compute_centre(array<shared_ptr<const Shell>, 2> shells) {

  const vector<double> exp0 = shells[0]->exponents();
  const vector<double> exp1 = shells[1]->exponents();
  array<double, 3> out = {{0.0, 0.0, 0.0}};
  for (auto& expi0 : exp0) {
    for (auto& expi1 : exp1) {
      const double cxp_inv = 1.0 / (expi0 + expi1);
      out[0] += (shells[0]->position(0) * expi0 + shells[1]->position(0) * expi1) * cxp_inv;
      out[1] += (shells[0]->position(1) * expi0 + shells[1]->position(1) * expi1) * cxp_inv;
      out[2] += (shells[0]->position(2) * expi0 + shells[1]->position(2) * expi1) * cxp_inv;
    }
  }
  const int denom = shells[0]->exponents().size() * shells[1]->exponents().size();
  out[0] /= denom;
  out[1] /= denom;
  out[2] /= denom;

  return out;
}


void Node::compute_multipoles(const int lmax) {

  const int nmultipole = (lmax + 1) * (lmax + 1);
  multipoles_.resize(nmultipole);
  vector<shared_ptr<ZMatrix>> multipoles(nmultipole);
  for (int i = 0; i != nmultipole; ++i)
    multipoles[i] = make_shared<ZMatrix>(nbasis_, nbasis_);

  bool skip = false;
  if (is_leaf_) { //compute multipole integrals
    size_t ob0 = 0;
    for (auto& a0 : bodies_) {
      for (auto& atom0 : a0->atoms()) {
        for (auto& b0 : atom0->shells()) {
          size_t ob1 = 0;
          for (auto& a1 : bodies_) {
            for (auto& atom1 : a1->atoms()) {
              for (auto& b1 : atom1->shells()) {
                //array<double, 3> centre = compute_centre({{b1, b0}});
                //MultipoleBatch mpole(array<shared_ptr<const Shell>, 2>{{b1, b0}}, centre, lmax);
                {
                  MultipoleBatch mpole(array<shared_ptr<const Shell>, 2>{{b1, b0}}, position_, lmax);
                  mpole.compute();
                  for (int i = 0; i != nmultipole; ++i)
                    multipoles[i]->copy_block(ob1, ob0, b1->nbasis(), b0->nbasis(), mpole.data(i));
                }

                ob1 += b1->nbasis();
              }
            }
          }
          ob0 += b0->nbasis();
        }
      }
    }
  } else { //shift children's multipoles
    shared_ptr<const Node> child = children_[0].lock();
    if (child->is_same_as_parent())
      skip = true;

    if (!skip) {
      size_t offset = 0;
      for (int n = 0; n != nchild_; ++n) {
        shared_ptr<const Node> child = children_[n].lock();
        array<double, 3> r12;
        r12[0] = position_[0] - child->position(0);
        r12[1] = position_[1] - child->position(1);
        r12[2] = position_[2] - child->position(2);
        assert(nmultipole == child->multipoles().size());
        LocalExpansion shift(r12, child->multipoles(), lmax);
        vector<shared_ptr<const ZMatrix>> moment = shift.compute_shifted_multipoles();

        for (int i = 0; i != nmultipole; ++i)
          multipoles[i]->copy_block(offset, offset, child->nbasis(), child->nbasis(), moment[i]->data());

        offset += child->nbasis();
      }
    }
  }

  if (!skip) {
    for (int i = 0; i != nmultipole; ++i)
      multipoles_[i] = multipoles[i];
  } else {
    shared_ptr<const Node> child = children_[0].lock();
    for (int i = 0; i != nmultipole; ++i)
      multipoles_[i] = child->multipoles(i);
  }
}


void Node::compute_local_expansions(shared_ptr<const Matrix> density, const int lmax, const vector<int> offsets) {

  const int nmultipole = (lmax + 1) * (lmax + 1);

  vector<int> l_map(nmultipole);
  int cnt = 0;
  for (int l = 0; l <= lmax; ++l)
    for (int m = 0; m <= 2 * l; ++m, ++cnt)
      l_map[cnt] = l;

  auto out = make_shared<ZMatrix>(nbasis_, nbasis_);

  if (density != nullptr) {
    for (auto& distant_node : interaction_list_) { // M2L
      array<double, 3> r12;
      r12[0] = position_[0] - distant_node->position(0);
      r12[1] = position_[1] - distant_node->position(1);
      r12[2] = position_[2] - distant_node->position(2);
      LocalExpansion lx(r12, distant_node->multipoles(), lmax);
      vector<shared_ptr<const ZMatrix>> lmoments = lx.compute_local_moments();

      // need to get a sub-density matrix corresponding to distant_node
      const int dimb = distant_node->nbasis();
      Matrix subden(dimb, dimb);
      subden.zero();
      size_t ob0 = 0;
      for (auto& body0 : distant_node->bodies()) {
        size_t iat0 = 0;
        for (auto& atom0 : body0->atoms()) {
          size_t ish0 = 0;
          for (auto& b0 : atom0->shells()) {
            const int offset0 = offsets[body0->ishell(iat0) + ish0];
            const size_t size0 = b0->nbasis();
            ++ish0;

            size_t ob1 = 0;
            for (auto& body1 : distant_node->bodies()) {
              size_t iat1 = 0;
              for (auto& atom1 : body1->atoms()) {
                size_t ish1 = 0;
                for (auto& b1 : atom1->shells()) {
                  const int offset1 = offsets[body1->ishell(iat1) + ish1];
                  const size_t size1 = b1->nbasis();
                  ++ish1;

                  shared_ptr<const Matrix> tmp = density->get_submatrix(offset1, offset0, size1, size0);
                  subden.copy_block(ob1, ob0, size1, size0, tmp);
                  ob1 += size1;
                }
                ++iat1;
              }
            }
            ob0 += size0;
          }
          ++iat0;
        }
      }

      for (int i = 0; i != nmultipole; ++i) {
        complex<double> contract = 0.0;
        for (int j = 0; j != dimb; ++j)
          for (int k = 0; k != dimb; ++k)
            contract += lmoments[i]->element(k, j) * subden.element(k, j);

        *out += pow(-1.0, l_map[i]) * contract * *multipoles_[i];
      }
    }
  }

  shared_ptr<const ZMatrix> nai = compute_NAI_far_field(lmax);
  local_expansion_ = make_shared<const ZMatrix>(*out + *nai);
  //local_expansion_ = make_shared<const ZMatrix>(*out);
}


shared_ptr<const ZMatrix> Node::compute_NAI_far_field(const int lmax) {

  const int nmultipole = (lmax + 1) * (lmax + 1);
  vector<int> l_map(nmultipole);
  int cnt = 0;
  for (int l = 0; l <= lmax; ++l)
    for (int m = 0; m <= 2 * l; ++m, ++cnt)
      l_map[cnt] = l;

  auto out = make_shared<ZMatrix>(nbasis_, nbasis_);

  for (auto& distant_node : interaction_list_) {
    array<double, 3> r12;
    r12[0] = distant_node->position(0) - position_[0];
    r12[1] = distant_node->position(1) - position_[1];
    r12[2] = distant_node->position(2) - position_[2];
    LocalExpansion lx(r12, multipoles_, lmax);
    vector<shared_ptr<const ZMatrix>> lmoments = lx.compute_local_moments();

    for (auto& body : distant_node->bodies()) {
      for (auto& atom : body->atoms()) {
        const double Z = atom->atom_charge();
        for (int i = 0; i != nmultipole; ++i)
          *out += -2.0 * pow(-1.0, l_map[i]) * Z * *lmoments[i];
      }
    }
  }

  return out;
}


shared_ptr<const ZMatrix> Node::compute_Coulomb(const int nbasis, shared_ptr<const Matrix> density, vector<int> offsets, const bool dodf, const string auxfile) {

  assert(is_leaf());
  auto out = make_shared<ZMatrix>(nbasis, nbasis);
  out->zero();

  // add FF local expansions to coulomb matrix
  size_t ob0 = 0;
  for (auto& body0 : bodies_) {
    size_t iat0 = 0;
    for (auto& atom0 : body0->atoms()) {
      size_t ish0 = 0;
      for (auto& b0 : atom0->shells()) {
        const int offset0 = offsets[body0->ishell(iat0) + ish0];
        const size_t size0 = b0->nbasis();
        ++ish0;

        size_t ob1 = 0;
        for (auto& body1 : bodies_) {
          size_t iat1 = 0;
          for (auto& atom1 : body1->atoms()) {
            size_t ish1 = 0;
            for (auto& b1 : atom1->shells()) {
              const int offset1 = offsets[body1->ishell(iat1) + ish1];
              const size_t size1 = b1->nbasis();
              ++ish1;

              ZMatrix sublocal = *(local_expansion_->get_submatrix(ob1, ob0, size1, size0));
              out->copy_block(offset1, offset0, size1, size0, sublocal.data());
              ob1 += size1;
            }
            ++iat1;
          }
        }
        ob0 += size0;
      }
      ++iat0;
    }
  }

  // compute near-field interactions using direct integration and add to far field
  vector<shared_ptr<const Shell>> basis;
  vector<int> new_offset;
  ////DEBUG
#if 0
  if (ninter_ != 0) {
    for (auto& body : bodies_) {
      size_t iat = 0;
      for (auto& atom : body->atoms()) {
        const vector<shared_ptr<const Shell>> tmp = atom->shells();
        basis.insert(basis.end(), tmp.begin(), tmp.end());
        vector<int> tmpoff;
        for (int ish = 0; ish != atom->shells().size(); ++ish)
          tmpoff.insert(tmpoff.end(), offsets[body->ishell(iat) + ish]);

        new_offset.insert(new_offset.end(), tmpoff.begin(), tmpoff.end());
        ++iat;
      }
    }
  }
#endif
  ////END OF DEBUG

  assert(iself_ >= 0);
  int ibas = -1;
  vector<shared_ptr<const Atom>> close_atoms;
  int nbas = 0;
  int inode = 0;
  for (auto& close_node : neighbour_) {
    if (inode == iself_) ibas = nbas;
    nbas += close_node->nbasis();
    for (auto& close_body : close_node->bodies()) {
      size_t iat = 0;
      for (auto& close_atom : close_body->atoms()) {
        close_atoms.push_back(close_atom);
        const vector<shared_ptr<const Shell>> tmp = close_atom->shells();
        basis.insert(basis.end(), tmp.begin(), tmp.end());
        vector<int> tmpoff;
        for (int ish = 0; ish != close_atom->shells().size(); ++ish)
          tmpoff.insert(tmpoff.end(), offsets[close_body->ishell(iat) + ish]);

        new_offset.insert(new_offset.end(), tmpoff.begin(), tmpoff.end());
        ++iat;
      }
    }
    ++inode;
  }
  assert (ibas >= 0);
  const size_t size = basis.size();

  // NAI for close-range
#if 1
  auto mol = make_shared<const Molecule>(close_atoms, vector<shared_ptr<const Atom>>{});

  for (auto& a3 : bodies_) {
    size_t iat3 = 0;
    for (auto& atom3 : a3->atoms()) {
      size_t ish3 = 0;
      for (auto& b3 : atom3->shells()) {
        const int b3size = b3->nbasis();
        const size_t b3offset = offsets[a3->ishell(iat3) + ish3];
        ++ish3;

        for (int i2 = 0; i2 != size; ++i2) {
          const shared_ptr<const Shell>  b2 = basis[i2];
          const int b2offset = new_offset[i2];
          const int b2size = b2->nbasis();
          {
            array<shared_ptr<const Shell>,2> shells = {{b3, b2}};
            NAIBatch naibatch(shells, mol);
            naibatch.compute();
            const double* naidata = naibatch.data();
            for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
              for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++naidata) {
                const double nai = *naidata;
                out->element(j2, j3) += 2.0 * nai;
              }
            }
          }
        }

      }
      ++iat3;
    }
  }
#endif


#if 1
  if (!dodf && density != nullptr) {
    const double* density_data = density->data();

    for (int i0 = 0; i0 != size; ++i0) {
      const shared_ptr<const Shell>  b0 = basis[i0];
      const int b0offset = new_offset[i0];
      const int b0size = b0->nbasis();

      for (int i1 = 0; i1 != size; ++i1) {
        const shared_ptr<const Shell>  b1 = basis[i1];
        const int b1offset = new_offset[i1];
        const int b1size = b1->nbasis();

        for (int i2 = 0; i2 != size; ++i2) {
          const shared_ptr<const Shell>  b2 = basis[i2];
          const int b2offset = new_offset[i2];
          const int b2size = b2->nbasis();

          for (auto& a3 : bodies_) {
            size_t iat3 = 0;
            for (auto& atom3 : a3->atoms()) {
              size_t ish3 = 0;
              for (auto& b3 : atom3->shells()) {
                const int b3size = b3->nbasis();
                const size_t b3offset = offsets[a3->ishell(iat3) + ish3];
                ++ish3;

                array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
                ERIBatch eribatch(input, 0.0);
                eribatch.compute();
                const double* eridata = eribatch.data();

                for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
                  const int j0n = j0 * density->ndim();
                  for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
                    for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                      for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                        const double eri = *eridata;
                        out->element(j2, j3) += density_data[j0n + j1] * eri;
                      }
                    }
                  }
                }

              }
              ++iat3;
            }
          }
        }
      }
    }
  } else if (dodf && density != nullptr) {
    vector<shared_ptr<const Atom>> aux_atoms;
    shared_ptr<const PTree> bdata = PTree::read_basis(auxfile);
    int naux =  0;
    for (auto& a : close_atoms) {
       auto aux_atom = make_shared<const Atom>(*a, a->spherical(), auxfile, make_pair(auxfile, bdata), nullptr);
       aux_atoms.push_back(aux_atom);
       naux += aux_atom->nbasis();
    }

    Matrix subden(nbas, nbas);
    int o0 = 0;
    for (int i0 = 0; i0 != size; ++i0) {
      const shared_ptr<const Shell>  b0 = basis[i0];
      const int b0offset = new_offset[i0];
      const int b0size = b0->nbasis();
      int o1 = 0;
      for (int i1 = 0; i1 != size; ++i1) {
        const shared_ptr<const Shell>  b1 = basis[i1];
        const int b1offset = new_offset[i1];
        const int b1size = b1->nbasis();
        auto tmp = density->get_submatrix(b1offset, b0offset, b1size, b0size);
        subden.copy_block(o1, o0, b1size, b0size, *tmp);

        o1 += b1size;
      }
      o0 += b0size;
    }

    shared_ptr<const DFDist> df = form_fit(nbas, naux, close_atoms, aux_atoms);
    shared_ptr<const Matrix> o = df->compute_Jop(make_shared<const Matrix>(subden));

    o0 = 0;
    for (int i0 = 0; i0 != size; ++i0) {
      const shared_ptr<const Shell>  b0 = basis[i0];
      const int b0offset = new_offset[i0];
      const int b0size = b0->nbasis();
      int o1 = ibas;
      for (auto& a1 : bodies_) {
        size_t iat1 = 0;
        for (auto& atom1 : a1->atoms()) {
          size_t ish1 = 0;
          for (auto& b1 : atom1->shells()) {
            const size_t b1offset = offsets[a1->ishell(iat1) + ish1];
            const int b1size = b1->nbasis();
            ++ish1;

            auto tmp = o->get_submatrix(o1, o0, b1size, b0size);
            out->add_real_block(1.0, b1offset, b0offset, b1size, b0size, *tmp);

            o1 += b1size;
          }
          ++iat1;
        }
      }
      o0 += b0size;
    }
  }
#endif

  return out;
}


shared_ptr<const ZMatrix> Node::compute_exact_Coulomb_FF(shared_ptr<const Matrix> density, vector<int> offsets) {

  assert(is_leaf());
  auto out = make_shared<ZMatrix>(density->ndim(), density->mdim());
  out->zero();

  if (ninter_ != 0) {
    vector<shared_ptr<const Shell>> basis;
    vector<int> new_offset;
    for (auto& body : bodies_) {
      size_t iat = 0;
      for (auto& atom : body->atoms()) {
        const vector<shared_ptr<const Shell>> tmp = atom->shells();
        basis.insert(basis.end(), tmp.begin(), tmp.end());
        vector<int> tmpoff;
        for (int ish = 0; ish != atom->shells().size(); ++ish)
          tmpoff.insert(tmpoff.end(), offsets[body->ishell(iat) + ish]);

        new_offset.insert(new_offset.end(), tmpoff.begin(), tmpoff.end());
        ++iat;
      }
    }
    for (auto& distant_node : interaction_list_) {
      for (auto& distant_body : distant_node->bodies()) {
        size_t iat = 0;
        for (auto& distant_atom : distant_body->atoms()) {
          const vector<shared_ptr<const Shell>> tmp = distant_atom->shells();
          basis.insert(basis.end(), tmp.begin(), tmp.end());
          vector<int> tmpoff;
          for (int ish = 0; ish != distant_atom->shells().size(); ++ish)
            tmpoff.insert(tmpoff.end(), offsets[distant_body->ishell(iat) + ish]);

          new_offset.insert(new_offset.end(), tmpoff.begin(), tmpoff.end());
          ++iat;
        }
      }
    }

    const size_t size = basis.size();
    const double* density_data = density->data();

    for (int i0 = 0; i0 != size; ++i0) {
      const shared_ptr<const Shell>  b0 = basis[i0];
      const int b0offset = new_offset[i0];
      const int b0size = b0->nbasis();

      for (int i1 = 0; i1 != size; ++i1) {
        const shared_ptr<const Shell>  b1 = basis[i1];
        const int b1offset = new_offset[i1];
        const int b1size = b1->nbasis();

        for (int i2 = 0; i2 != size; ++i2) {
          const shared_ptr<const Shell>  b2 = basis[i2];
          const int b2offset = new_offset[i2];
          const int b2size = b2->nbasis();

          for (auto& a3 : bodies_) {
            size_t iat3 = 0;
            for (auto& atom3 : a3->atoms()) {
              size_t ish3 = 0;
              for (auto& b3 : atom3->shells()) {
                const int b3size = b3->nbasis();
                const size_t b3offset = offsets[a3->ishell(iat3) + ish3];
                ++ish3;

                array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
                ERIBatch eribatch(input, 0.0);
                eribatch.compute();
                const double* eridata = eribatch.data();

                for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
                  const int j0n = j0 * density->ndim();
                  for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
                    for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                      for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                        const double eri = *eridata;
                        out->element(j2, j3) += density_data[j0n + j1] * eri;
                      }
                    }
                  }
                }
              }
              ++iat3;
            }
          }
        }
      }
    }
  }

  return out;
}


void Node::print_node() const {

  cout << "*** Node at " << setprecision(5) << position_[0] << "   " << position_[1] << "   " << position_[2] << endl;
  cout << "has " << nbody_ << " bodies: " << endl;
  for (auto& body : bodies_)
    cout << setw(5) << setprecision(5) << setw(10) << body->position(0) << "   "
                                       << setw(10) << body->position(1) << "   "
                                       << setw(10) << body->position(2) << endl;
  cout << "and " << nneighbour_ << " neighbours" << endl;
  for (auto& neigh : neighbour_)
    cout << setprecision(5) << setw(10) << neigh->position(0) << "  "
                            << setw(10) << neigh->position(1) << "  "
                            << setw(10) << neigh->position(2) << endl;
}
