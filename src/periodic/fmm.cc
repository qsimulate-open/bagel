//
// BAGEL - Parallel electron correlation program.
// Filename: fmm.cc
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


#include <src/periodic/fmm.h>
#include <src/util/taskqueue.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/periodic/multipolebatch.h>
#include <src/periodic/localexpansion.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;
const static Legendre plm;

FMM::FMM(shared_ptr<const Geometry> geom, const int ns, const int lmax, const double thresh, const int ws)
 : geom_(geom), ns_(ns), lmax_(lmax), thresh_(thresh), ws_(ws) {

  init();
}


void FMM::init() {

  centre_ = geom_->charge_center();
  nbasis_ = geom_->nbasis();
  const int ns2 = pow(2, ns_);

  nsp_ = geom_->nshellpair();

  double maxext = 0;
  for (int i = 1; i != nsp_; ++i)
    if (maxext < geom_->shellpair(i)->extent())
      maxext = geom_->shellpair(i)->extent();
  const int maxws = maxext / ws_;
  if (maxws > ns2) throw runtime_error("maxws > 2**ns");

  coordinates_.resize(nsp_);
  maxxyz_ = {{0, 0, 0}};
  double rad = 0;
  for (int j = 0; j != 3; ++j) {
    for (int i = 0; i != nsp_; ++i) {
      coordinates_[i][j] = geom_->shellpair(i)->centre(j);
      if (abs(coordinates_[i][j]) > maxxyz_[j]) maxxyz_[j] = abs(coordinates_[i][j]);
    }
    if (maxxyz_[j] > rad) rad = maxxyz_[j];
  }
  boxsize_  = 2.05 * rad;
  unitsize_ = boxsize_/ns2;

  cout << "boxsize = " << boxsize_ << " unitsize = " << unitsize_ << " maxxyz = " << maxxyz_[0] << " " << maxxyz_[1] << " " << maxxyz_[2] << endl;

  get_boxes();
}


void FMM::get_boxes() {

  Timer fmminit;

  const int ns2 = pow(2, ns_);

  // find out unempty leaves
  vector<array<int, 3>> boxid; // max dimension 2**(ns+1)-1
  vector<int> ibox(nsp_);

  map<array<int, 3>, int> treemap;
  assert(treemap.empty());
  int nleaf = 0;
  for (int isp = 0; isp != nsp_; ++isp) {
    array<int, 3> idxbox;
    for (int i = 0; i != 3; ++i) {
      idxbox[i] = (int) floor(coordinates_[isp][i]/unitsize_) + ns2/2;
      assert(idxbox[i] < ns2);
    }

    map<array<int, 3>,int>::iterator box = treemap.find(idxbox);
    const bool box_found = (box != treemap.end());
    if (box_found) {
      ibox[isp] = treemap.find(idxbox)->second;
    } else {
      treemap.insert(treemap.end(), pair<array<int, 3>,int>(idxbox, nleaf));
      ibox[isp] = nleaf;
      boxid.resize(nleaf+1);
      boxid[nleaf] = idxbox;
      ++nleaf;
    }
  }
  assert(nleaf == boxid.size() && nleaf <= nsp_);

  // get leaves
  vector<vector<int>> leaves(nleaf);
  for (int isp = 0; isp != nsp_; ++isp) {
    const int n = ibox[isp];
    assert(n < nleaf);
    leaves[n].insert(leaves[n].end(), isp);
  }

  // get all unempty boxes
  int nbox = 0;
  for (int il = 0; il != nleaf; ++il) {
    vector<shared_ptr<const ShellPair>> sp;
    for (int i = 0; i != leaves[il].size(); ++i)
      sp.insert(sp.end(), geom_->shellpair(leaves[il][i]));
    auto newbox = make_shared<Box>(0, il, boxid[il], lmax_, sp);
    box_.insert(box_.end(), newbox);
    ++nbox;
  }

  int icntc = 0;
  int icntp = ns2;
  nbranch_.resize(ns_+2);
  nbranch_[0] = nleaf;

  for (int nss = ns_; nss > -1; --nss) {
    int nbranch = 0;
    const int nss2 = pow(2, nss);

    for (int i = 0; i != nss2; ++i) {
      for (int j = 0; j != nss2; ++j) {
        for (int k = 0; k != nss2; ++k) {
          vector<shared_ptr<const ShellPair>> sp;
          array<int, 3> idxp;
          idxp[0] = (int) floor(0.5*(i+1)) + icntp;
          idxp[1] = (int) floor(0.5*(j+1)) + icntp;
          idxp[2] = (int) floor(0.5*(k+1)) + icntp;

          array<int, 3> idxc = {{i+icntc, j+icntc, k+icntc}};
          map<array<int, 3>,int>::iterator child = treemap.find(idxc);
          const bool child_found = (child != treemap.end());

          if (child_found) {
            const int ichild = treemap.find(idxc)->second;
            map<array<int, 3>,int>::iterator parent = treemap.find(idxp);
            const bool parent_found = (parent != treemap.end());

            if (!parent_found) {
              auto newbox = make_shared<Box>(ns_-nss+1, nbox, idxp, lmax_, box_[ichild]->sp());
              box_.insert(box_.end(), newbox);
              treemap.insert(treemap.end(), pair<array<int, 3>,int>(idxp, nbox));
              box_[nbox]->insert_child(box_[ichild]);
              ++nbox;
              ++nbranch;
            } else {
              const int ibox = treemap.find(idxp)->second;
              box_[ibox]->insert_child(box_[ichild]);
              box_[ibox]->insert_sp(box_[ichild]->sp());
            }
          }

        }
      }
    }
    icntc = icntp;
    icntp += nss2;
    nbranch_[ns_-nss+1] = nbranch;
  }
  assert(accumulate(nbranch_.begin(), nbranch_.end(), 0) == nbox);
  nbox_ = nbox;
  cout << "ns_ = " << ns_ << " nbox = " << nbox_ << "  nleaf = " << nleaf << " nsp = " << nsp_ << endl;

  for (auto& b : box_)
    b->init();

  int icnt = 0;
  for (int ir = 0; ir != ns_+1; ++ir) {
    vector<shared_ptr<Box>> tmpbox(nbranch_[ir]);
    for (int ib = 0; ib != nbranch_[ir]; ++ib)
      tmpbox[ib] = box_[icnt + ib];
    for (auto& b : tmpbox)
      b->get_neigh(tmpbox, ws_);
    icnt += nbranch_[ir];
  }

  fmminit.tick_print("fmm initialisation");

}


void FMM::M2M() {

  Timer m2mtime;
  m2mtime.tick_print("Upward M2M pass");
}


void FMM::M2L() {

  Timer dtime;
  dtime.tick_print("M2L pass");
}


void FMM::L2L() {

  Timer dtime;
  dtime.tick_print("Downward L2L pass");
}


shared_ptr<const ZMatrix> FMM::compute_energy(shared_ptr<const Matrix> density, const bool dodf, const double schwarz_thresh) const {

  const size_t nbasis = density->ndim();
  auto out = make_shared<ZMatrix>(nbasis, nbasis);

  return out;
}


void FMM::print_boxes(const int i) const {

  int ib = 0;
  for (auto& b : box_) {
    if (b->rank() == i) {
      cout << "Box " << ib << " rank = " << i << " *** nchild = " << b->nchild() << " *** nsp = " << b->nsp() << " *** Shell pairs at:" << endl;
      for (int i = 0; i != b->nsp(); ++i)
        cout << setprecision(5) << b->sp(i)->centre(0) << "  " << b->sp(i)->centre(1) << "  " << b->sp(i)->centre(2) << endl;
      ++ib;
    }
    if (b->rank() > i) break;
  }

}
