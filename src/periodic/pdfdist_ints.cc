//
// BAGEL - Parallel electron correlation program.
// Filename: pdfdist_ints.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/periodic/pdfdist_ints.h>
#include <src/periodic/pdfinttask.h>

using namespace bagel;
using namespace std;

PDFDist_ints::PDFDist_ints(const vector<array<double, 3>>& L, const int nbas, const int naux,
                           const vector<shared_ptr<const Atom>>& atoms_c0,
                           const vector<shared_ptr<const Atom>>& atoms_cg,
                           const vector<shared_ptr<const Atom>>& aux_atoms, const shared_ptr<const Geometry> cell0,
                           const double thr, const shared_ptr<const Matrix> p, const shared_ptr<const VectorB> data1)
  : DFDist(nbas, naux), lattice_vectors_(L), projector_(p), data1_(data1) {

  assert(nindex1_ == nindex2_);
  nbasis_ = nindex1_;

  // 3index integrals made in DFBlock.
  vector<shared_ptr<const Shell>> ashell, b0shell, bgshell;
  for (auto& i : aux_atoms)     ashell.insert(ashell.end(),  i->shells().begin(), i->shells().end());
  for (auto& i : atoms_c0)     b0shell.insert(b0shell.end(), i->shells().begin(), i->shells().end());
  for (auto& i : atoms_cg)     bgshell.insert(bgshell.end(), i->shells().begin(), i->shells().end());

  int astart;
  vector<shared_ptr<const Shell>> myashell;
  tie(astart, myashell) = get_ashell(ashell);

  shared_ptr<const StaticDist> adist_shell = make_table(astart);
  shared_ptr<const StaticDist> adist_averaged = make_shared<const StaticDist>(naux_, mpi__->size());

  // make an empty dfblock for (rs|i)
  const size_t asize  = accumulate(myashell.begin(),myashell.end(),0, [](const int& i, const shared_ptr<const Shell>& o) { return i+o->nbasis(); });
  const size_t b0size = accumulate(b0shell.begin(), b0shell.end(), 0, [](const int& i, const shared_ptr<const Shell>& o) { return i+o->nbasis(); });
  const size_t bgsize = accumulate(bgshell.begin(), bgshell.end(), 0, [](const int& i, const shared_ptr<const Shell>& o) { return i+o->nbasis(); });
  assert(naux_ == asize && b0size == nbasis_ && bgsize == nbasis_);
  block_.push_back(make_shared<DFBlock>(adist_shell, adist_averaged, asize, b0size, bgsize, astart, 0, 0));
  block_[0]->zero();

  // 3-index integrals (r sL'|iL) for each L' (sum over all lattice vectors L)
  pcompute_3index(myashell, b0shell, bgshell);

  // compute NAI
  pcompute_NAI(b0shell, bgshell, cell0);

  // charged part of coeff
  compute_charged_coeff(b0shell, bgshell);
}


void PDFDist_ints::pcompute_3index(const vector<shared_ptr<const Shell>>& ashell,
                                   const vector<shared_ptr<const Shell>>& b0shell,
                                   const vector<shared_ptr<const Shell>>& bgshell) {
  Timer time;

  TaskQueue<PDFIntTask_3index> tasks(b0shell.size() * bgshell.size() * ashell.size() * ncell());
  auto i3 = make_shared<const Shell>(ashell.front()->spherical());

  data3_in_cell_.resize(ncell());
  for (auto& block : data3_in_cell_)
    block = make_shared<DFBlock>(*block_[0]);

  int j2 = 0;
  for (auto& i2 : b0shell) {
    int j1 = 0;
    for (auto& i1 : bgshell) {
      int n = 0;
      for (auto& L : lattice_vectors_) {
        int j0 = 0;
        for (auto& i0 : ashell) {
          auto i00 = make_shared<const Shell>(*(i0->move_atom(L)));
          tasks.emplace_back((array<shared_ptr<const Shell>, 4>{{i3, i00, i1, i2}}), (array<int, 3>{{j2, j1, j0}}), data3_in_cell_[n]);
          j0 += i0->nbasis();
        }
        ++n;
      }
      j1 += i1->nbasis();
    }
    j2 += i2->nbasis();
  }

  tasks.compute();
//time.tick_print("3-index integrals");

  // now project and sum
  if (!projector_)
    throw logic_error("failed attempt to project data3_ without the projection matrix");

  // P_C = 1 - P
  auto projectorC = make_shared<Matrix>(naux_, naux_, serial_);
  projectorC->unit();
  *projectorC -= *projector_;

  auto tmp = make_shared<btas::Tensor3<double>>(naux_, nbasis_, nbasis_);
  assert(block_[0]->size() == naux_ * nbasis_ * nbasis_);
  for (int i = 0; i != ncell(); ++i) {
    contract(1.0, *data3_in_cell_[i], {3, 1, 2}, *projectorC, {0, 3}, 0.0, *tmp, {0, 1, 2});
    blas::ax_plus_y_n(1.0, tmp->data(), block_[0]->size(), block_[0]->data());
  }

}


void PDFDist_ints::pcompute_NAI(const vector<shared_ptr<const Shell>>& b0shell,
                                const vector<shared_ptr<const Shell>>& bgshell,
                                const shared_ptr<const Geometry> cell0) {
  Timer time;

  TaskQueue<PDFIntTask_NAI> tasks(b0shell.size() * bgshell.size() * ncell());

  nai_in_cell_.resize(ncell());
  for (auto& mat : nai_in_cell_) mat = make_shared<Matrix>(nbasis_, nbasis_);

  int j1 = 0;
  for (auto& i1 : bgshell) {
    int j0 = 0;
    for (auto& i0 : b0shell) {
      int n = 0;
      for (auto& L : lattice_vectors_) {
        auto mol = make_shared<const Geometry>(*cell0, L);
        tasks.emplace_back((array<shared_ptr<const Shell>, 2>{{i1, i0}}), (array<int, 2>{{j0, j1}}), mol, nai_in_cell_[n]);
        ++n;
      }
      j0 += i0->nbasis();
    }
    j1 += i1->nbasis();
  }

  tasks.compute();
//time.tick_print("NAI integrals");

}


void PDFDist_ints::compute_charged_coeff(const vector<shared_ptr<const Shell>>& b0shell,
                                         const vector<shared_ptr<const Shell>>& bgshell) {

  Timer time;
  coeffC_ = make_shared<btas::Tensor3<double>>(naux_, nbasis_, nbasis_);
  fill_n(coeffC_->data(), naux_*nbasis_*nbasis_, 0.0);

  if (!data1_)
    throw logic_error("auxiliary charge has to be computed first");

  // first compute <t|uL>
  TaskQueue<PDFIntTask_coeff> tasks(b0shell.size() * bgshell.size());

  int j1 = 0;
  for (auto& b1 : bgshell) {
    int j0 = 0;
    for (auto& b0 : b0shell) {
      tasks.emplace_back((array<shared_ptr<const Shell>, 2>{{b1, b0}}), (array<int, 2>{{j0, j1}}), coeffC_, data1_);
      j0 += b0->nbasis();
    }
    j1 += b1->nbasis();
  }

  tasks.compute();
//time.tick_print("overlap integrals");
}
