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


#include <src/util/f77.h>
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
  nshell_ = 0;
  double sum = 0.0;
  for (auto& body : bodies_) {
    for (auto& atom : body->atoms()) {
      position_[0] += atom->atom_charge() * atom->position(0);
      position_[1] += atom->atom_charge() * atom->position(1);
      position_[2] += atom->atom_charge() * atom->position(2);
      sum += atom->atom_charge();
      nbasis_ += atom->nbasis();
      nshell_ += atom->nshell();
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


bool Node::is_neighbour(array<shared_ptr<const Shell>, 4> shells, const int ws) {

  const double extent01 = pair_extent(array<shared_ptr<const Shell>, 2>{{shells[0], shells[1]}});
  const double extent23 = pair_extent(array<shared_ptr<const Shell>, 2>{{shells[2], shells[3]}});

  array<double, 3> rvec01 = compute_centre(array<shared_ptr<const Shell>, 2>{{shells[0], shells[1]}});;
  array<double, 3> rvec23 = compute_centre(array<shared_ptr<const Shell>, 2>{{shells[2], shells[3]}});;

  array<double, 3> v;
  v[0] = rvec01[0] - rvec23[0];
  v[1] = rvec01[1] - rvec23[1];
  v[2] = rvec01[2] - rvec23[2];
  const double r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return (r <= (1.0 + ws) * (extent01 + extent23));
}


double Node::pair_extent(array<shared_ptr<const Shell>, 2> shells) {

  const vector<double> exp0 = shells[0]->exponents();
  const vector<double> exp1 = shells[1]->exponents();
  array<double, 3> AB;
  AB[0] = shells[0]->position(0) - shells[1]->position(0);
  AB[1] = shells[0]->position(1) - shells[1]->position(1);
  AB[2] = shells[0]->position(2) - shells[1]->position(2);
  const double rsq = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
  const double lnthresh = log(thresh_);

  double out = 0;
  array<double, 3> rvec = compute_centre(array<shared_ptr<const Shell>, 2>{{shells[0], shells[1]}});;
  for (auto& expi0 : exp0) {
    for (auto& expi1 : exp1) {
      const double cxp_inv = 1.0 / (expi0 + expi1);
      const double expi01 = expi0 * expi1;
      const double lda_kl = sqrt(abs(- lnthresh - expi01 * rsq * cxp_inv + 0.75 * log(4.0 * expi01 / pisq__)) * cxp_inv);
      //const double s01 = pow(4.0 * expi01 * cxp_inv * cxp_inv, 0.75) * exp(-expi01 * cxp_inv * rsq);
      //const double r01 = sqrt((-lnthresh + log(s01) + 0.5 * log(expi0 + expi1)) * cxp_inv);

      array<double, 3> tmp;
      tmp[0] = (shells[0]->position(0) * expi0 + shells[1]->position(0) * expi1) * cxp_inv - rvec[0];
      tmp[1] = (shells[0]->position(1) * expi0 + shells[1]->position(1) * expi1) * cxp_inv - rvec[1];
      tmp[2] = (shells[0]->position(2) * expi0 + shells[1]->position(2) * expi1) * cxp_inv - rvec[2];

      const double extent01 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]) + lda_kl;
      //const double extent01 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]) + r01;
      if (extent01 > out) out = extent01;
    }
  }

  return out;
}


void Node::compute_extent(const double thresh) {

  vector<shared_ptr<const Shell>> shells;
  for (auto& body : bodies_)
    for (auto& atom : body->atoms())
      shells.insert(shells.end(), atom->shells().begin(), atom->shells().end());

  extent_ = 0.0;
  for (auto& ish : shells)
    for (auto& jsh : shells) {
      const double extent01 = pair_extent(array<shared_ptr<const Shell>, 2>{{ish, jsh}});
      if (extent01 > extent_) extent_ = extent01;
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
        sort_neighbours(neighbour_);
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
    sort_neighbours(neighbour_);
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


void Node::compute_local_expansions(shared_ptr<const Matrix> density, const int lmax, const vector<int> offsets, const double scale) {

  const int nmultipole = (lmax + 1) * (lmax + 1);

  vector<int> l_map(nmultipole);
  int cnt = 0;
  for (int l = 0; l <= lmax; ++l)
    for (int m = 0; m <= 2 * l; ++m, ++cnt)
      l_map[cnt] = l;

  int nfbas = 0;
  for (auto& distant_node : interaction_list_)
    nfbas += distant_node->nbasis();

  auto lrs = make_shared<ZMatrix>(nbasis_, nbasis_);
  auto lrt = make_shared<ZMatrix>(nbasis_, nfbas);

  if (density) {
    size_t ob = 0;
    for (auto& distant_node : interaction_list_) { // M2L
      array<double, 3> r12;
      r12[0] = position_[0] - distant_node->position(0);
      r12[1] = position_[1] - distant_node->position(1);
      r12[2] = position_[2] - distant_node->position(2);
      LocalExpansion lx(r12, distant_node->multipoles(), lmax);
      vector<shared_ptr<const ZMatrix>> lmoments = lx.compute_local_moments();

      // get sub-density matrix D_tu (tu=far)
      const int dimb = distant_node->nbasis();
      Matrix den_tu(dimb, dimb);
      den_tu.zero();
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
                  den_tu.copy_block(ob1, ob0, size1, size0, tmp);
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
        const double contract = lmoments[i]->get_real_part()->dot_product(den_tu);
        *lrs += pow(-1.0, l_map[i]) * contract * *multipoles_[i];
      }

#if 0
      // get sub-density matrix D_su (s=this u=far)
      auto den_su = make_shared<Matrix>(nbasis_, dimb);
      ob0 = 0;
      for (auto& body0 : distant_node->bodies()) {
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

                  shared_ptr<const Matrix> tmp = density->get_submatrix(offset1, offset0, size1, size0);
                  den_su->copy_block(ob1, ob0, size1, size0, tmp);
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

      auto zden_su = make_shared<ZMatrix>(*den_su, 1.0);
      auto tmp_ts = make_shared<ZMatrix>(nbasis_, dimb);
      auto tmp_rt = make_shared<ZMatrix>(nbasis_, dimb);
      for (int i = 0; i != nmultipole; ++i) {
        zgemm3m_("N", "N", nbasis_, dimb, dimb, 1.0, zden_su->data(), nbasis_, lmoments[i]->data(), dimb, 0.0, tmp_ts->data(), nbasis_);
        zgemm3m_("N", "N", nbasis_, dimb, nbasis_, 1.0, multipoles_[i]->data(), nbasis_, tmp_ts->data(), nbasis_, 0.0, tmp_rt->data(), nbasis_);
        const int sign = pow(-1.0, l_map[i]);
        lrt->add_block(sign, 0, ob, nbasis_, dimb, tmp_rt);
      }
#endif

      ob += distant_node->nbasis();
    }
  }

//  shared_ptr<const ZMatrix> nai = compute_NAI_far_field(lmax, scale);
//  local_expansion_ = make_shared<const ZMatrix>(*lrs + *nai);
  local_expansion_ = make_shared<const ZMatrix>(*lrs);
  //exchange_ = make_shared<const ZMatrix>(*lrt);
}


shared_ptr<const ZMatrix> Node::compute_NAI_far_field(const int lmax, const double scale) {

  auto out = make_shared<ZMatrix>(nbasis_, nbasis_);

  for (auto& distant_node : interaction_list_) {
    array<double, 3> r12;
    r12[0] = distant_node->position(0) - position_[0];
    r12[1] = distant_node->position(1) - position_[1];
    r12[2] = distant_node->position(2) - position_[2];
    LocalExpansion lx(r12, multipoles_, lmax);
    vector<shared_ptr<const ZMatrix>> lmoments = lx.compute_local_moments();

    for (auto& body : distant_node->bodies())
      for (auto& atom : body->atoms())
        *out += -scale * atom->atom_charge() * *lmoments[0];
  }

  return out;
}


shared_ptr<const ZMatrix> Node::compute_Coulomb(const int nbasis, shared_ptr<const Matrix> density, vector<int> offsets, const bool dodf, const double scale, const vector<double> schwarz, const double schwarz_thresh) {

  assert(is_leaf());
  auto out = make_shared<ZMatrix>(nbasis, nbasis);
  out->zero();
  n1e_int_ = 0;
  n2e_int_ = 0;
  n2e_total_ = 0;
  const size_t ndim = density->ndim();

#if 0
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

#if 0
        size_t ob = 0;
        for (auto& distant_node : interaction_list_) {
          for (auto& body1 : distant_node->bodies()) {
            size_t iat1 = 0;
            for (auto& atom1 : body1->atoms()) {
              size_t ish1 = 0;
              for (auto& b1 : atom1->shells()) {
                const int offset1 = offsets[body1->ishell(iat1) + ish1];
                const size_t size1 = b1->nbasis();
                ++ish1;

                ZMatrix sublocal = *(exchange_->get_submatrix(ob0, ob, size0, size1));
                out->copy_block(offset0, offset1, size0, size1, sublocal.data());
                ob += size1;
              }
              ++iat1;
            }
          }
        }
#endif

        ob0 += size0;
      }
      ++iat0;
    }
  }
#endif


#if 1 //DEBUG: no near-field
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
        for (int ish = 0; ish != atom->shells().size(); ++ish)
          new_offset.insert(new_offset.end(), offsets[body->ishell(iat) + ish]);
        ++iat;
      }
    }
  }
#endif
  ////END OF DEBUG

  assert(iself_ >= 0);
  int ibas = -1;
  int ishell = -1;
  vector<shared_ptr<const Atom>> close_atoms;
  vector<int> shell_id;
  int nbas = 0;
  int nsh = 0;
  int inode = 0;
  for (auto& close_node : neighbour_) {
    if (inode == iself_) {
      ibas = nbas;
      ishell = nsh;
    }
    nbas += close_node->nbasis();
    for (auto& close_body : close_node->bodies()) {
      size_t iat = 0;
      for (auto& close_atom : close_body->atoms()) {
        close_atoms.push_back(close_atom);
        const vector<shared_ptr<const Shell>> tmp = close_atom->shells();
        basis.insert(basis.end(), tmp.begin(), tmp.end());
        for (int ish = 0; ish != tmp.size(); ++ish) {
          new_offset.insert(new_offset.end(), offsets[close_body->ishell(iat) + ish]);
          shell_id.insert(shell_id.end(), close_body->ishell(iat) + ish);
        }

        nsh += close_atom->nshell();
        ++iat;
      }
    }
    ++inode;
  }
  assert(ibas >= 0 && ishell >= 0);
  const size_t size = basis.size();
  assert(ishell < size);
  assert(shell_id.size() == size);


#if 0 /////// DEBUG
  // get significant distributions
  int nsignif = 0;
  for (int i = 0; i != size; ++i) {
    const shared_ptr<const Shell>  bi = basis[i];
    const int ioffset = new_offset[i];
    const int isize = bi->nbasis();
    for (int j = 0; j != size; ++j) {
      const shared_ptr<const Shell>  bj = basis[j];
      const int joffset = new_offset[j];
      const int jsize = bj->nbasis();

      {
        array<shared_ptr<const Shell>,2> shells = {{bi, bj}};
        OverlapBatch obatch(shells);
        obatch.compute();
        const double* odata = obatch.data();
        double tmp = 0;
        for (int n = 0; n != isize * jsize; ++n, ++odata)
          tmp = *odata * *odata;
        if (tmp > 1e-10) ++nsignif;
      }
    }
  }
//  cout << size * size << " *** " << nsignif << endl;
#endif ///////DEBUG

#if 0
  // NAI for close-range
  auto mol = make_shared<const Molecule>(close_atoms, vector<shared_ptr<const Atom>>{});
  n1e_int_ = nshell_ * size;

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
            array<shared_ptr<const Shell>,2> shells = {{b2, b3}};
            NAIBatch naibatch(shells, mol);
            naibatch.compute();
            const double* naidata = naibatch.data();
            for (int j3 = b3offset; j3 != b3offset + b3size; ++j3) {
              for (int j2 = b2offset; j2 != b2offset + b2size; ++j2, ++naidata) {
                const double nai = *naidata;
                out->element(j3, j2) += scale * nai;
              }
            }
          }
        }

      }
      ++iat3;
    }
  }
#endif

  n2e_total_ = nshell_ * size * size * size;
  if (!dodf && density) {
    shellquads_.clear();
    int_offsets_.clear();
    const double* density_data = density->data();
    const int nshell = sqrt(schwarz.size());

    vector<double> max_density(size * size);
    vector<int> do_int(size * size, 1);
    for (int i = 0; i != size; ++i) {
      const int ioffset = new_offset[i];
      const int isize = basis[i]->nbasis();
      for (int j = i; j != size; ++j) {
        const int joffset = new_offset[j];
        const int jsize = basis[j]->nbasis();

        const int ij = i * size + j;
        const int ji = j * size + i;

        const double schwarz_ij = schwarz[shell_id[i]*nshell + shell_id[j]];
        if (schwarz_ij < 1e-15) {
          do_int[ij] = 0;
          do_int[ji] = 0;
        }


        double denmax = 0.0;
        for (int ii = ioffset; ii != ioffset + isize; ++ii) {
          const int iin = ii * density->ndim();
          for (int jj = joffset; jj != joffset + jsize; ++jj)
            denmax = max(denmax, fabs(density_data[iin + jj]));
        }
        max_density[ij] = denmax;
        max_density[ji] = denmax;
      }
    }

    int nfar = 0;
    for (int i0 = 0; i0 != size; ++i0) {
      const shared_ptr<const Shell>  b0 = basis[i0];
      const size_t b0offset = new_offset[i0];
      const size_t b0size = b0->nbasis();

      for (int i1 = 0; i1 != size; ++i1) {
//        if (shell_id[i1] < shell_id[i0]) continue;
        const shared_ptr<const Shell>  b1 = basis[i1];
        const size_t b1offset = new_offset[i1];
        const size_t b1size = b1->nbasis();

        const int i01 = i0 * size + i1;
        const double density_01 = max_density[i01] * 4.0;
        const double schwarz_i01 = schwarz[shell_id[i0]*nshell + shell_id[i1]];

        const bool skip01 = do_int[i0 * size + i1] == 0;
        if (skip01) continue;

        for (int i2 = 0; i2 != size; ++i2) {
//          if (shell_id[i2] < shell_id[i0]) continue;
          const shared_ptr<const Shell>  b2 = basis[i2];
          const size_t b2offset = new_offset[i2];
          const size_t b2size = b2->nbasis();

          const double density_02 = max_density[i0 * size + i2];
          const double density_12 = max_density[i1 * size + i2];

          int i3 = ishell;
          for (auto& a3 : bodies_) {
            size_t iat3 = 0;
            for (auto& atom3 : a3->atoms()) {
              size_t ish3 = 0;
              for (auto& b3 : atom3->shells()) {
                const size_t b3size = b3->nbasis();
                const size_t b3offset = offsets[a3->ishell(iat3) + ish3];
                ++ish3;
//                if (shell_id[i3] < shell_id[i2]) continue;

                const bool skip23 = do_int[i2 * size + i3] == 0;
                if (skip23) continue;

                const bool isclose = is_neighbour(array<shared_ptr<const Shell>, 4>{{b0, b1, b2, b3}}, 2);
                if (!isclose) {
                  ++nfar;
                  continue;
                }

                const int i23 = i2 * size + i3;
//                if (shell_id[i2]*nshell + shell_id[i3] < shell_id[i0]*nshell + shell_id[i1]) continue;
                const double density_23 = max_density[i23] * 4.0;
                const double density_03 = max_density[i0 * size + i3];
                const double density_13 = max_density[i1 * size + i3];
                const double schwarz_i23 = schwarz[shell_id[i2]*nshell + shell_id[i3]];

                const double mulfactor = max(max(max(density_01, density_02),
                                                 max(density_12, density_23)),
                                                 max(density_03, density_13));
                const double integral_bound = mulfactor * schwarz_i01 * schwarz_i23;
                const bool skip_schwarz = integral_bound < schwarz_thresh;
                if (skip_schwarz) continue;
                ++n2e_int_;

                shellquads_.insert(shellquads_.end(), {{shell_id[i0], shell_id[i1], shell_id[i2], shell_id[i3]}});
                int_offsets_.insert(int_offsets_.end(), {{b0offset, b1offset, b2offset, b3offset}});
#if 0
/////////// DEBUG //////// NO NF INTEGRATION
                array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
                ERIBatch eribatch(input, mulfactor);
                eribatch.compute();
                const double* eridata = eribatch.data();

                for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
                  const int j0n = j0 * density->ndim();
                  for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
                    const int j1n = j1 * density->ndim();
                    for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                      for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                        double eri = *eridata;
                        out->element(j2, j3) += density_data[j0n + j1] * eri;
                        out->element(j0, j3) -= density_data[j1n + j2] * eri * 0.5;
                      }
                    }
                  }
                }
#endif

                ++i3;
              }
              ++iat3;
            }
          }
        }
      }
    }
  } else if (dodf && density) {
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

    shared_ptr<const Matrix> o = df_->compute_Jop(make_shared<const Matrix>(subden));

#if 1
    // exchange
    shared_ptr<Matrix> coeff = subden.copy();
    *coeff *= -1.0;
    int nocc = 0;
    {
      VectorB vec(subden.ndim());
      coeff->diagonalize(vec);
      for (int i = 0; i != subden.ndim(); ++i) {
        if (vec[i] < -1.0e-8) {
          ++nocc;
          const double fac = std::sqrt(-vec(i));
          for_each(coeff->element_ptr(0,i), coeff->element_ptr(0,i+1), [&fac](double& i) { i *= fac; });
        } else { break; }
      }
    }
    if (nocc == 0) return out;
    shared_ptr<DFHalfDist> halfbj = df_->compute_half_transform(coeff->slice(0,nocc));
    shared_ptr<DFHalfDist> half = halfbj->apply_J();
    Matrix ex = *half->form_2index(half, -0.5);
#endif

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

            auto tmp = o->get_submatrix(o0, o1, b0size, b1size);
            out->add_real_block(1.0, b0offset, b1offset, b0size, b1size, *tmp); // J

            auto tmpex = ex.get_submatrix(o0, o1, b0size, b1size);
            out->add_real_block(1.0, b0offset, b1offset, b0size, b1size, *tmpex); // K

            o1 += b1size;
          }
          ++iat1;
        }
      }
      o0 += b0size;
    }
  }
#endif //DEBUG: no near-field

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
        for (int ish = 0; ish != atom->shells().size(); ++ish)
          new_offset.insert(new_offset.end(), offsets[body->ishell(iat) + ish]);
        ++iat;
      }
    }
    for (auto& distant_node : interaction_list_) {
      for (auto& distant_body : distant_node->bodies()) {
        size_t iat = 0;
        for (auto& distant_atom : distant_body->atoms()) {
          const vector<shared_ptr<const Shell>> tmp = distant_atom->shells();
          basis.insert(basis.end(), tmp.begin(), tmp.end());
          for (int ish = 0; ish != distant_atom->shells().size(); ++ish)
            new_offset.insert(new_offset.end(), offsets[distant_body->ishell(iat) + ish]);
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


void Node::sort_neighbours(vector<shared_ptr<const Node>> neighbours) {

  const size_t nneigh = neighbours.size();
  const size_t nbit = neighbours[0]->key().size();
  vector<shared_ptr<const Node>>::iterator it = neighbours.begin();
  for (int i = 0; i != nneigh - 1; ++i) {
    bool swap = false;
    for (int j = nbit-1; j >= 0; j--)
      if (neighbours[i]->key()[j] ^ neighbours[nneigh-1]->key()[j])
        swap = neighbours[nneigh-1]->key()[j];
    if (swap) {
      neighbours.insert(it, neighbours[nneigh-1]);
      neighbours.resize(nneigh);
      break;
    }
    ++it;
  }
}


void Node::form_df(const string auxfile) {

  vector<shared_ptr<const Atom>> close_atoms;
  int nbas = 0;
  for (auto& close_node : neighbour_) {
    nbas += close_node->nbasis();
    for (auto& close_body : close_node->bodies())
      for (auto& close_atom : close_body->atoms())
      close_atoms.insert(close_atoms.end(), close_atom);
  }

  vector<shared_ptr<const Atom>> aux_atoms;
  shared_ptr<const PTree> bdata = PTree::read_basis(auxfile);
  int naux =  0;
  for (auto& a : close_atoms) {
     auto aux_atom = make_shared<const Atom>(*a, a->spherical(), auxfile, make_pair(auxfile, bdata), nullptr);
     aux_atoms.push_back(aux_atom);
     naux += aux_atom->nbasis();
  }
  df_ = form_fit(nbas, naux, close_atoms, aux_atoms);
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
