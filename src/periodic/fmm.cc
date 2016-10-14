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
#include <src/periodic/jexpansion.h>
#include <boost/math/special_functions/erf.hpp>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;
static const double tolgd__ = 1e-9;
const static Legendre plm;

FMM::FMM(shared_ptr<const Geometry> geom, const int ns, const int lmax, const double thresh, const int ws)
 : geom_(geom), ns_(ns), lmax_(lmax), thresh_(thresh), ws_(ws) {

  init();
}


void FMM::init() {

  centre_ = geom_->charge_center();
  cout << "**** CHARGE CENTRE " << setprecision(9) << centre_[0] << "  " << centre_[1] << "   " << centre_[2] << endl;
  nbasis_ = geom_->nbasis();
  const int ns2 = pow(2, ns_);

  int nsh = 0;
  for (auto& a : geom_->atoms())
    nsh += a->shells().size();
  cout << "******** NSHELL = " << nsh << " *********** " << endl;

  const int nsp = geom_->nshellpair();
  base_extent_ = sqrt(2.0) * boost::math::erf_inv(1.0-tolgd__);

#if 0
  double maxext = 0;
  for (int i = 1; i != nsp_; ++i)
    if (maxext < geom_->shellpair(i)->extent())
      maxext = geom_->shellpair(i)->extent();
  const int maxws = maxext / ws_;
  if (maxws > ns2) throw runtime_error("maxws > 2**ns");
#endif

  maxxyz_ = {{0, 0, 0}};
  double rad = 0;
  int isp = 0;
  for (int i = 0; i != nsp; ++i) {
    if (geom_->shellpair(i)->extent() < numerical_zero__) continue;
    isp_.push_back(i);
    ++isp;
    for (int j = 0; j != 3; ++j) {
      const double jco = geom_->shellpair(i)->centre(j);
      if (abs(jco) > maxxyz_[j]) maxxyz_[j] = abs(jco);
      if (maxxyz_[j] > rad) rad = maxxyz_[j];
    }
  }
  assert(isp_.size() == isp);
  nsp_ = isp;
  coordinates_.resize(nsp_);
  for (int i = 0; i != nsp_; ++i)
    coordinates_[i] = geom_->shellpair(isp_[i])->centre();

  boxsize_  = 2.05 * rad;
  unitsize_ = boxsize_/ns2;
  coordinates_.resize(nsp_);
  isp_.resize(nsp_);

  cout << "boxsize = " << boxsize_ << " unitsize = " << unitsize_ << " maxxyz = " << maxxyz_[0] << " " << maxxyz_[1] << " " << maxxyz_[2] << endl;

  get_boxes();

  do_ff_ = false;
  for (int i = 0; i != nbranch_[0]; ++i)
    if (box_[i]->ninter() != 0) do_ff_ = true;
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
      const double coi = coordinates_[isp][i]-centre_[i];
      idxbox[i] = (int) floor(coi/unitsize_) + ns2/2 + 1;
      assert(idxbox[i] <= ns2 && idxbox[i] > 0);
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
    const int ingeom = isp_[isp];
    leaves[n].insert(leaves[n].end(), ingeom);
  }

  // get all unempty boxes
  int nbox = 0;
  for (int il = 0; il != nleaf; ++il) {
    vector<shared_ptr<const ShellPair>> sp;
    for (auto& isp : leaves[il])
      sp.insert(sp.end(), geom_->shellpair(isp));
    array<int, 3> id = boxid[il];
    auto newbox = make_shared<Box>(0, il, id, lmax_, sp);
    box_.insert(box_.end(), newbox);
    ++nbox;
  }

  int icntc = 0;
  int icntp = ns2;
  nbranch_.resize(ns_+1);

  for (int nss = ns_; nss > -1; --nss) {
    int nbranch = 0;
    const int nss2 = pow(2, nss);

    for (int i = 1; i != nss2+1; ++i) {
      for (int j = 1; j != nss2+1; ++j) {
        for (int k = 1; k != nss2+1; ++k) {
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
              if (nss != 0) {
                auto newbox = make_shared<Box>(ns_-nss+1, nbox, idxp, lmax_, box_[ichild]->sp());
                box_.insert(box_.end(), newbox);
                treemap.insert(treemap.end(), pair<array<int, 3>,int>(idxp, nbox));
                box_[nbox]->insert_child(box_[ichild]);
                box_[ichild]->insert_parent(box_[nbox]);
                ++nbox;
              }
              ++nbranch;
            } else {
              const int ibox = treemap.find(idxp)->second;
              box_[ibox]->insert_child(box_[ichild]);
              box_[ibox]->insert_sp(box_[ichild]->sp());
              box_[ichild]->insert_parent(box_[ibox]);
              ++nbranch;
            }
          }

        }
      }
    }
    icntc = icntp;
    icntp += nss2;
    nbranch_[ns_-nss] = nbranch;
  }
  assert(accumulate(nbranch_.begin(), nbranch_.end(), 0) == nbox);
  nbox_ = nbox;
  cout << "ns_ = " << ns_ << " nbox = " << nbox_ << "  nleaf = " << nleaf << " nsp = " << nsp_ << endl;

  for (auto& b : box_)
    b->init();

  int icnt = 0;
  for (int ir = ns_; ir > -1; --ir) {
    vector<shared_ptr<Box>> tmpbox(nbranch_[ir]);
    for (int ib = 0; ib != nbranch_[ir]; ++ib)
      tmpbox[ib] = box_[nbox_-icnt-nbranch_[ir]+ib];
    for (auto& b : tmpbox) {
      b->get_neigh(tmpbox, ws_);
      b->get_inter(tmpbox, ws_);
    }
    icnt += nbranch_[ir];
  }

  int i = 0;
  for (auto& b : box_) {
    const bool ipar = (b->parent()) ? true : false;
    cout << i << " rank = " << b->rank() << " extent = " << b->extent() << " nsp = " << b->nsp()
         << " nchild = " << b->nchild() << " nneigh = " << b->nneigh() << " ninter = " << b->ninter()
         << " centre = " << b->centre(0) << " " << b->centre(1) << " " << b->centre(2)
         << " idxc = " << b->tvec()[0] << " " << b->tvec()[1] << " " << b->tvec()[2] << " *** " << ipar << endl;
    ++i;
  }

  fmminit.tick_print("fmm initialisation");
}


void FMM::M2M(shared_ptr<const Matrix> density) const {

  Timer m2mtime;
  int u = 0;
  for (int i = 0; i != nbranch_[0]; ++i)
    if (u++ % mpi__->size() == mpi__->rank())
      box_[i]->compute_M2M(density);
  m2mtime.tick_print("shift sp");

  //u = 0;
  int icnt = nbranch_[0];
  for (int i = 1; i != ns_+1; ++i)
  //  if (u++ % mpi__->size() == mpi__->rank())
      for (int j = 0; j != nbranch_[i]; ++j, ++icnt)
        box_[icnt]->compute_M2M(density);

  assert(icnt == nbox_);

  m2mtime.tick_print("M2M pass");
}


void FMM::M2L() const {

  Timer m2ltime;

  for (auto& b : box_)
    b->compute_M2L();

  m2ltime.tick_print("M2L pass");
}


void FMM::L2L() const {

  Timer l2ltime;

  int icnt = 0;
  for (int ir = ns_; ir > -1; --ir) {
    for (int ib = 0; ib != nbranch_[ir]; ++ib)
      box_[nbox_-icnt-nbranch_[ir]+ib]->compute_L2L();

    icnt += nbranch_[ir];
  }

  l2ltime.tick_print("L2L pass");
}


shared_ptr<const ZMatrix> FMM::compute_energy(shared_ptr<const Matrix> density) const {

  Timer nftime;
  auto out = make_shared<ZMatrix>(nbasis_, nbasis_);
  out->zero();
 
  M2M(density);
  M2L();
  L2L();

  if (density) {
    assert(nbasis_ == density->ndim());
    vector<double> maxden(geom_->nshellpair());
    const double* density_data = density->data();
    for (int i01 = 0; i01 != geom_->nshellpair(); ++i01) {
      shared_ptr<const Shell> sh0 = geom_->shellpair(i01)->shell(1);
      const int offset0 = geom_->shellpair(i01)->offset(1);
      const int size0 = sh0->nbasis();

      shared_ptr<const Shell> sh1 = geom_->shellpair(i01)->shell(0);
      const int offset1 = geom_->shellpair(i01)->offset(0);
      const int size1 = sh0->nbasis();

      double denmax = 0.0;
      for (int i0 = offset0; i0 != offset0 + size0; ++i0) {
        const int i0n = i0 * density->ndim();
        for (int i1 = offset1; i1 != offset1 + size1; ++i1)
          denmax = max(denmax, fabs(density_data[i0n + i1]));
      }
      maxden[i01] = denmax;
    }

    int u = 0;
    for (int i = 0; i != nbranch_[0]; ++i) {
      if (u++ % mpi__->size() == mpi__->rank()) {
        auto ei = box_[i]->compute_node_energy(density, maxden, geom_->schwarz_thresh());
        *out += *ei;
      }
    }
    out->allreduce();

    for (int i = 0; i != nbasis_; ++i) out->element(i, i) *= 2.0;
    out->fill_upper();
  }

  nftime.tick_print("near-field");

  return out;
}


double FMM::energy_ff() const {

  double out = 0;
  for (int i = 0; i != nbranch_[0]; ++i) {
    box_[i]->compute_energy_ff();
    out += box_[i]->energy_ff();
  }

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
