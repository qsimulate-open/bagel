//
// BAGEL - Parallel electron correlation program.
// Filename: node_sp.cc
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


#include <src/util/f77.h>
#include <src/periodic/node_sp.h>
#include <src/periodic/multipolebatch.h>
#include <src/integral/os/overlapbatch.h>
#include <src/periodic/localexpansion.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;

NodeSP::NodeSP(const bitset<nbit__> key, const int depth, shared_ptr<const NodeSP> parent, const double thresh)
 : key_(key), depth_(depth), parent_(parent), thresh_(thresh) {
  if (depth == 0)
    key_[0] = 1;

  is_leaf_ = false;
  is_complete_ = false;
  nchild_  = 0;
  nvertex_   = 0;
  nneigh_ = 0;
  ninter_ = 0;
  is_same_as_parent_ = false;
  rank_ = 0;
}


void NodeSP::insert_vertex(shared_ptr<const VertexSP> p) {

  assert(!is_complete_);
  vertex_.resize(nvertex_ + 1);
  vertex_[nvertex_] = p;

  ++nvertex_;
}


void NodeSP::insert_child(shared_ptr<const NodeSP> child) {

  if (child) {
    assert(!is_leaf_);
    child_.resize(nchild_ + 1);
    child_[nchild_] = child;
    ++nchild_;
  }
}


void NodeSP::init() {

  is_complete_ = true;
  extent_ = -1.0;
  if (nchild_ == 0) is_leaf_ = true;

  if (!is_leaf_)
    child_local_expansion_.resize(nchild_);

  nbasis0_ = 0;
  nbasis1_ = 0;
  centre_ = {{0.0, 0.0, 0.0}};
  for (auto& v : vertex_) {
    centre_[0] += v->centre(0);
    centre_[1] += v->centre(1);
    centre_[2] += v->centre(2);
    nbasis0_ += v->nbasis0();
    nbasis1_ += v->nbasis1();
  }
  centre_[0] /= nvertex_;
  centre_[1] /= nvertex_;
  centre_[2] /= nvertex_;

  array<double, 3> r12;
  r12[0] = centre_[0] - parent_->centre(0);
  r12[1] = centre_[1] - parent_->centre(1);
  r12[2] = centre_[2] - parent_->centre(2);

  const double r = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
  if (r <= numerical_zero__) {
    is_same_as_parent_ = true;
    rank_ = parent_->rank();
  } else {
    rank_ = parent_->rank() + 1;
  }

  compute_extent();
}

double NodeSP::radius() const {

  double out = 0;
  for (auto& v : vertex_) {
    array<double, 3> r12;
    r12[0] = centre_[0] - v->centre(0);
    r12[1] = centre_[1] - v->centre(1);
    r12[2] = centre_[2] - v->centre(2);
    const double r = sqrt(r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2]);
    if (r > out) out = r;
  }

  return out;
}

bool NodeSP::is_neighbour(shared_ptr<const NodeSP> nd, const int ws) {

  array<double, 3> rvec = nd->centre();
  assert(extent_ > 0 && nd->extent() > 0);

  array<double, 3> v;
  v[0] = centre_[0] - rvec[0];
  v[1] = centre_[1] - rvec[1];
  v[2] = centre_[2] - rvec[2];
  const double r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return (r <= (1.0 + ws) * (extent_ + nd->extent()));
}


void NodeSP::compute_extent() {

  extent_ = 0.0;
  vector<shared_ptr<const Shell>> shells;
  for (auto& v : vertex_)
    if (v->extent() > extent_) extent_ = v->extent();
}


void NodeSP::insert_neigh(shared_ptr<const NodeSP> neigh, const bool is_neigh, const int ws) {

  assert(neigh->depth() == depth_);
  if (!is_neigh) {
    const bool check = is_neighbour(neigh, ws);
    if (check) {
      neigh_.resize(nneigh_ + 1);
      neigh_[nneigh_] = neigh;
      ++nneigh_;
    } else {
      interaction_list_.resize(ninter_ + 1);
      interaction_list_[ninter_] = neigh;
      ++ninter_;
    }
  } else {
    neigh_.resize(nneigh_ + 1);
    neigh_[nneigh_] = neigh;
    ++nneigh_;
  }
}


void NodeSP::make_interaction_list(const int ws) {

  ninter_ = 0;
  for (auto& pneigh : parent_->neigh()) {
    for (int n = 0; n != pneigh->nchild(); ++n) {
      shared_ptr<const NodeSP> child = pneigh->child(n);
      assert(child->depth() == depth_);
      bool is_close = false;
      for (int i = 0; i != nneigh_; ++i) {
        if (child->key() == neigh_[i]->key()) {
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


void NodeSP::compute_multipoles(const int lmax) {

  const int nmult = (lmax + 1) * (lmax + 1);
  multipole_.resize(nmult);
  vector<shared_ptr<ZMatrix>> multipole(nmult);
  for (int i = 0; i != nmult; ++i)
    multipole[i] = make_shared<ZMatrix>(nbasis0_, nbasis1_);

  bool skip = false;
  if (is_leaf_) { // shift sp's multipoles
    size_t offset0 = 0;
    size_t offset1 = 0;
    for (auto& v : vertex_) {
      array<double, 3> r12;
      r12[0] = centre_[0] - v->centre(0);
      r12[1] = centre_[1] - v->centre(1);
      r12[2] = centre_[2] - v->centre(2);
      assert(nmult == v->nmult());
      const double r = sqrt(r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2]);
      //cout << setprecision(9) << " r = " << r << endl;
//      LocalExpansion shift(r12, v->multipole(), lmax);
//      vector<shared_ptr<const ZMatrix>> moment = shift.compute_shifted_multipoles();

      for (int i = 0; i != nmult; ++i)
        //multipole[i]->copy_block(offset1, offset0, v->nbasis1(), v->nbasis0(), moment[i]->data());
        multipole[i]->copy_block(offset1, offset0, v->nbasis1(), v->nbasis0(), v->multipole()[i]->data());

      offset0 += v->nbasis0();
      offset1 += v->nbasis1();
    }
    assert(offset0 == nbasis0_);
    assert(offset1 == nbasis1_);
  } else { // shift children's multipoles
    shared_ptr<const NodeSP> child = child_[0].lock();
    if (child->is_same_as_parent())
      skip = true;

    if (!skip) {
      size_t offset0 = 0;
      size_t offset1 = 0;
      for (int n = 0; n != nchild_; ++n) {
        shared_ptr<const NodeSP> child = child_[n].lock();
        array<double, 3> r12;
        r12[0] = centre_[0] - child->centre(0);
        r12[1] = centre_[1] - child->centre(1);
        r12[2] = centre_[2] - child->centre(2);
        assert(nmult == child->multipole().size());
        LocalExpansion shift(r12, child->multipole(), lmax);
        vector<shared_ptr<const ZMatrix>> moment = shift.compute_shifted_multipoles();

        for (int i = 0; i != nmult; ++i)
          multipole[i]->copy_block(offset1, offset0, child->nbasis1(), child->nbasis0(), moment[i]->data());

        offset0 += child->nbasis0();
        offset1 += child->nbasis1();
      }
      assert(offset0 == nbasis0_);
      assert(offset1 == nbasis1_);
    }
  }

  if (!skip) {
    for (int i = 0; i != nmult; ++i)
      multipole_[i] = multipole[i];
  } else {
    shared_ptr<const NodeSP> child = child_[0].lock();
    for (int i = 0; i != nmult; ++i)
      multipole_[i] = child->multipole(i);
  }
}


void NodeSP::compute_local_expansions(shared_ptr<const Matrix> density, const int lmax) {

  if (lmax < 0 || !density) return;

  const int nmult = (lmax + 1) * (lmax + 1);
  vector<int> l_map(nmult);
  int cnt = 0;
  for (int l = 0; l <= lmax; ++l)
    for (int m = 0; m <= 2 * l; ++m, ++cnt)
      l_map[cnt] = l;

  auto lrs = make_shared<ZMatrix>(nbasis1_, nbasis0_);

  if (density) {
    for (auto& distant_node : interaction_list_) { // M2L
      array<double, 3> r12;
      r12[0] = centre_[0] - distant_node->centre(0);
      r12[1] = centre_[1] - distant_node->centre(1);
      r12[2] = centre_[2] - distant_node->centre(2);
      LocalExpansion lx(r12, distant_node->multipole(), lmax);
      vector<shared_ptr<const ZMatrix>> lmoments = lx.compute_local_moments();

      // get sub-density matrix D_tu (tu=far)
      const int dimb0 = distant_node->nbasis0();
      const int dimb1 = distant_node->nbasis1();
      assert(lmoments[0]->ndim() == dimb1 && lmoments[0]->mdim() == dimb0);
      Matrix den_tu(dimb1, dimb0);
      den_tu.zero();
      size_t ob0 = 0;
      size_t ob1 = 0;
      for (auto& v : distant_node->vertex()) {

        shared_ptr<const Shell> sh0 = v->shell0();
        shared_ptr<const Shell> sh1 = v->shell1();

        const int size0 = sh0->nbasis();
        const int size1 = sh1->nbasis();

        shared_ptr<const Matrix> den = density->get_submatrix(v->offset(1), v->offset(0), size1, size0);
        den_tu.copy_block(ob1, ob0, size1, size0, den);

#if 0
        /******** DEBUG: contract with multipoles here and compare with 4c integrals ********/
        if (is_leaf_) {

          // get mult of v
          vector<shared_ptr<const ZMatrix>> vmult(nmult);
          for (int i = 0; i != nmult; ++i)
            vmult[i] = distant_node->multipole(i)->get_submatrix(ob1, ob0, size1, size0);

          // check mult of v
          vector<shared_ptr<const ZMatrix>> vchk(nmult);
          {
            MultipoleBatch mpole(v->shells(), v->centre(), lmax);
            mpole.compute();

            for (int i = 0; i != nmult; ++i) {
              ZMatrix tmp0(size1, size0);
              tmp0.copy_block(0, 0, size1, size0, mpole.data(i));
              vchk[i] = make_shared<const ZMatrix>(tmp0);
            }
          }
#if 1
          for (int i = 0; i != nmult; ++i) {
            auto verr = make_shared<const ZMatrix>(*vchk[i] - *vmult[i]);
            const double vrmserr = verr->get_real_part()->rms();
            if (vrmserr > 1e-2) {
              cout << " *** CHECK MULT OF V *** L = " << i << " error = " << setprecision(9) << vrmserr << endl;
              vmult[i]->get_real_part()->print("MULT FROM DISTANT NODE");
              vchk[i]->get_real_part()->print("MULT FROM BATCH");
            }
          }
#endif


          // get local expansion of v (using r12 to distant node)
          LocalExpansion local(r12, vchk, lmax);
          //LocalExpansion local(r12, vmult, lmax);
          vector<shared_ptr<const ZMatrix>> vl = local.compute_local_moments();

          // for each v0, compute jrs using mult of v0
          size_t o0 = 0;
          size_t o1 = 0;
          vector<shared_ptr<const ZMatrix>> v0mult(nmult);
          for (auto& v0 : vertex_) {
            const double r = sqrt(r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2]);
            // check if v and v0 are well sep
            const bool isneigh = v->is_neighbour(v0, 2);

            shared_ptr<const Shell> s0 = v0->shell0();
            shared_ptr<const Shell> s1 = v0->shell1();
            for (int i = 0; i != nmult; ++i)
              v0mult[i] = multipole_[i]->get_submatrix(o1, o0, s1->nbasis(), s0->nbasis());

            // check mult of v0
            vector<shared_ptr<const ZMatrix>> vchk0(nmult);
            {
              MultipoleBatch mpole0(v0->shells(), v0->centre(), lmax);
              mpole0.compute();

              for (int i = 0; i != nmult; ++i) {
                ZMatrix tmp0(s1->nbasis(), s0->nbasis());
                tmp0.copy_block(0, 0, s1->nbasis(), s0->nbasis(), mpole0.data(i));
                vchk0[i] = make_shared<const ZMatrix>(tmp0);
              }
            }

            auto jrs = make_shared<ZMatrix>(s1->nbasis(), s0->nbasis());
            for (int i = 0; i != nmult; ++i) {
              const double contract = vl[i]->get_real_part()->dot_product(den);
              *jrs += pow(-1.0, l_map[i]) * contract * *vchk0[i];
              //*jrs += pow(-1.0, l_map[i]) * contract * *v0mult[i];
            }

            // now compare jrs with 4c-jrs using v and v0
            auto jrs0 = make_shared<ZMatrix>(s1->nbasis(), s0->nbasis());
            array<shared_ptr<const Shell>,4> input = {{sh1, sh0, s1, s0}};

            {
              ERIBatch eribatch(input, 0.0);
              eribatch.compute();
              const double* eridata = eribatch.data();

              for (int j0 = 0; j0 != s0->nbasis(); ++j0)
                for (int j1 = 0; j1 != s1->nbasis(); ++j1)
                  for (int j2 = 0; j2 != size0; ++j2)
                    for (int j3 = 0; j3 != size1; ++j3, ++eridata) {
                      const double eri = *eridata;
                      jrs0->element(j1, j0) += den->element(j2, j3) * eri;
                    }
            }

#if 0
            auto error = make_shared<const ZMatrix>(*jrs-*jrs0);
            const double rmserr = error->get_real_part()->rms();
            if (rmserr > 1e-4) {
              cout << " **CHECK MULT APPROX *** " << setprecision(9) << error->get_real_part()->rms() << "  R = " << r << " check neigh = " << isneigh << endl;
              jrs0->get_real_part()->print("4C INTEGRATION");
              jrs->get_real_part()->print("MULTIPOLE APPROXIMATION");
            }
#endif

            o0 += s0->nbasis();
            o1 += s1->nbasis();
          }
        }
        /******** END DEBUG ********/
#endif

        ob1 += size1;
        ob0 += size0;
      }
      assert(ob1 == dimb1);
      assert(ob0 == dimb0);

      for (int i = 0; i != nmult; ++i) {
        const double contract = lmoments[i]->get_real_part()->dot_product(den_tu);
        *lrs += pow(-1.0, l_map[i]) * contract * *multipole_[i];
      }

    }
  }

  local_expansion_ = make_shared<const ZMatrix>(*lrs);
}


shared_ptr<const ZMatrix> NodeSP::compute_Coulomb(const int dim, shared_ptr<const Matrix> density, vector<double> max_den, const bool dodf, const double schwarz_thresh) {

  assert(is_leaf());
  auto out = make_shared<ZMatrix>(dim, dim);
  out->zero();
  const size_t ndim = density->ndim();

  // add FF local expansions to coulomb matrix
  #if 0
  if (local_expansion_) {
    size_t ob0 = 0;
    size_t ob1 = 0;
    for (auto& v : vertex_) {
      shared_ptr<const Shell> sh0 = v->shell0();
      shared_ptr<const Shell> sh1 = v->shell1();
      const int size0 = sh0->nbasis();
      const int size1 = sh1->nbasis();

      ZMatrix sublocal = *(local_expansion_->get_submatrix(ob1, ob0, size1, size0));
      out->copy_block(v->offset(1), v->offset(0), size1, size0, sublocal.data());

      ob0 += size0;
      ob1 += size1;
    }
  }
  #endif

  // compute near-field interactions using direct integration and add to far field

  if (!dodf && density) {
    const double* density_data = density->data();
    const int nsh = sqrt(max_den.size());

    for (auto& v01 : vertex_) {
//      if (v01->schwarz() < 1e-15) continue;
      shared_ptr<const Shell> b0 = v01->shell0();
      const int i0 = v01->shell_ind(0);
      const int b0offset = v01->offset(0);
      const int b0size = b0->nbasis();

      shared_ptr<const Shell> b1 = v01->shell1();
      const int i1 = v01->shell_ind(1);
      const int b1offset = v01->offset(1);
      const int b1size = b1->nbasis();

      const int i01 = i0 * density->ndim() + i1;
      const double density_01 = max_den[i01] * 4.0;

      for (auto& neigh : neigh_) {
//        if (neigh->id_in_tree() < id_in_tree_) continue;
        for (auto& v23 : neigh->vertex()) {
//          if (v23->schwarz() < 1e-15) continue;
          shared_ptr<const Shell> b2 = v23->shell0();
          const int i2 = v23->shell_ind(0);
          const int b2offset = v23->offset(0);
          const int b2size = b2->nbasis();

          shared_ptr<const Shell> b3 = v23->shell1();
          const int i3 = v23->shell_ind(1);
          const int b3offset = v23->offset(1);
          const int b3size = b3->nbasis();

          const int i23 = i2 * density->ndim() + i3;

          const double density_23 = max_den[i01] * 4.0;
          const double density_02 = max_den[i0 * nsh + i2];
          const double density_03 = max_den[i0 * nsh + i3];
          const double density_12 = max_den[i1 * nsh + i2];
          const double density_13 = max_den[i1 * nsh + i3];

          const double mulfactor = max(max(max(density_01, density_02),
                                           max(density_03, density_12)),
                                           max(density_13, density_23));
          const double integral_bound = mulfactor * v01->schwarz() * v23->schwarz();
          const bool skip_schwarz = integral_bound < schwarz_thresh;
          //if (skip_schwarz) continue;

          array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
          ERIBatch eribatch(input, mulfactor);
          eribatch.compute();
          const double* eridata = eribatch.data();

          for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
            const int j0n = j0 * density->ndim();
            for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
              const int j1n = j1 * density->ndim();
              for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                const int j2n = j2 * density->ndim();
                for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                  const int j3n = j3 * density->ndim();
                  //const double scale = ((neigh->id_in_tree() == id_in_tree_) ? 1.0 : 2.0);

                  double eri = *eridata;
                  out->element(j0, j1) += density_data[j2n + j3] * eri;
                  out->element(j0, j2) -= density_data[j1n + j3] * eri * 0.5;
                  //if (neigh->id_in_tree() !=  id_in_tree_) {
                  //  out->element(j2, j3) += density_data[j0n + j1] * eri;
                  //  out->element(j1, j3) -= density_data[j0n + j2] * eri * 0.5;
                  //}
                }
              }
            }
          }

        }
      }
    }
  } else if (dodf && density) {
    throw runtime_error("No DF yet :-(");
  }

  return out;
}


void NodeSP::sort_neigh() {

  vector<shared_ptr<const NodeSP>> tmp(neigh_.size());

  const size_t nbit = neigh_[0]->key().size();
  vector<shared_ptr<const NodeSP>>::iterator it = neigh_.begin();
  for (int i = 0; i != nneigh_-1; ++i) {
    bool swap = false;
    for (int j = nbit-1; j >= 0; j--)
      if (neigh_[i]->key()[j] ^ neigh_[nneigh_-1]->key()[j])
        swap = neigh_[nneigh_-1]->key()[j];

    if (swap) {
      neigh_.insert(it, neigh_[nneigh_-1]);
      neigh_.resize(nneigh_);
      break;
    }
    ++it;
  }
}


void NodeSP::print_node() const {

  cout << "*** Node at " << setprecision(5) << centre_[0] << "   " << centre_[1] << "   " << centre_[2] << endl;
  cout << "has " << nvertex_ << " vertex: " << endl;
  for (auto& vertex : vertex_)
    cout << setw(5) << setprecision(5) << setw(10) << vertex->centre(0) << "   "
                                       << setw(10) << vertex->centre(1) << "   "
                                       << setw(10) << vertex->centre(2) << endl;
  cout << "and " << nneigh_ << " neighs" << endl;
  for (auto& neigh : neigh_)
    cout << setprecision(5) << setw(10) << neigh->centre(0) << "  "
                            << setw(10) << neigh->centre(1) << "  "
                            << setw(10) << neigh->centre(2) << endl;
}
