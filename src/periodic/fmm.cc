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

  maxxyz_ = {{0, 0, 0}};
  double rad = 0;
  for (int j = 0; j != 3; ++j) {
    for (int i = 0; i != nsp_; ++i) {
      coordinates_[i][j] = geom_->shellpair(i)->centre(j) - centre_[j];
      if (coordinates_[i][j] > maxxyz_[j]) maxxyz_[j] = coordinates_[i][j];
    }
    if (maxxyz_[j] > rad) rad = maxxyz_[j];
  }
  boxsize_  = 2.05 * rad;
  unitsize_ = boxsize_/ns2;
  nsp_box_.resize(nsp_);

  get_boxes();
}


void FMM::get_boxes() {

  const int ns2 = pow(2, ns_);

  // find out unempty leaves
  vector<array<int, 3>> boxid; // max dimension 2**(ns+1)-1
  vector<int> ibox(nsp_);

  map<array<int, 3>, int> leafmap;
  int nleaf = 0;
  for (int isp = 0; isp != nsp_; ++isp) {
    array<int, 3> idxbox;
    for (int i = 0; i != 3; ++i) {
      idxbox[i] = (int) floor(coordinates_[isp][i]/unitsize_);
      assert(idxbox[i] < ns2);
    }

    pair<map<array<int, 3>,int>::iterator, bool> chk_newbox;
    chk_newbox = leafmap.insert(pair<array<int, 3>,int>(idxbox, nleaf));
    if (!chk_newbox.second) {
     leafmap.insert(leafmap.end(), pair<array<int, 3>,int>(idxbox, nleaf));
     ibox[isp] = nleaf;
     boxid[nleaf] = idxbox;
     ++nleaf;
    } else {
      ibox[isp] = leafmap.find(idxbox)->second;
    }
  }
  assert(nleaf == boxid.size());

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
      sp.insert(sp.end(), geom_->shellpair(i));
    auto newbox = make_shared<Box>(il, boxid[il], lmax_, sp);
    box_.insert(box_.end(), newbox);
    ++nbox;
  }

  int icntc = 0;
  int icntp = ns2;
  for (int nss = ns_; nss != 0; --nss) {
    const int nss2 = pow(nss, 2);

    for (int i = 0; i != nss2; ++i) {
      for (int j = 0; i != nss2; ++j) {
        for (int k = 0; i != nss2; ++k) {
          vector<shared_ptr<const ShellPair>> sp;
          array<int, 3> idxp;
          idxp[0] = (int) floor(0.5*(i+1)) + icntp;
          idxp[1] = (int) floor(0.5*(j+1)) + icntp;
          idxp[2] = (int) floor(0.5*(k+1)) + icntp;
          array<int, 3> idxc = {{i+icntc, j+icntc, k+icntc}};

          pair<map<array<int, 3>,int>::iterator, bool> child_exist;
          child_exist = leafmap.insert(pair<array<int, 3>,int>(idxc, 0));
          pair<map<array<int, 3>,int>::iterator, bool> parent_exist;
          parent_exist = leafmap.insert(pair<array<int, 3>,int>(idxp, 0));
          if (child_exist.second && !parent_exist.second) {
            const int ichild = leafmap.find(idxc)->second;
            auto newbox = make_shared<Box>(nbox, idxp, lmax_, box_[ichild]->sp());
            box_.insert(box_.end(), newbox);
            leafmap.insert(leafmap.end(), pair<array<int, 3>,int>(idxp, nbox));
            ++nbox;
          } else if (child_exist.second && parent_exist.second) {
            const int ibox = leafmap.find(idxp)->second;
            const int ichild = leafmap.find(idxc)->second;
            box_[ibox]->insert_sp(box_[ichild]->sp());
          }

        }
      }
    }
    icntc = icntp;
    icntp += nss2;
  }

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
