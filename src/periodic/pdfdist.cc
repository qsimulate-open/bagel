//
// BAGEL - Parallel electron correlation program.
// Filename: pdfdist.cc
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

#include <src/periodic/pdfdist.h>
#include <src/periodic/pdfinttask.h>

using namespace bagel;
using namespace std;

PDFDist::PDFDist(vector<array<double, 3>> L, const int nbas, const int naux,
                 const vector<shared_ptr<const Atom>>& atoms0,
                 const vector<shared_ptr<const Atom>>& aux_atoms,
                 const double thresh, const bool serial, const shared_ptr<Matrix> data2)
  : lattice_vectors_(L), nbasis_(nbas), naux_(naux), serial_(serial), data2_(data2) {

  // prepare to compute 1- and 2-index ints
  vector<shared_ptr<const Shell>> ashell, b0shell;
  for (auto& i : atoms0)   b0shell.insert(b0shell.end(), i->shells().begin(), i->shells().end());
  for (auto& i : aux_atoms) ashell.insert(ashell.end(),  i->shells().begin(), i->shells().end());

  // compute auxiliary charge data1_ then projection matrix
  compute_aux_charge(ashell);

  // 2-index integrals (i|j_L)^{-1} (sum over L)
  if (data2)
    data2_ = data2;
  else
    pcompute_2index(ashell, thresh);

  /** form object PDFDist_ints for every cell, each cell L contains the 3-index integral
      (r0sL|aL') sum over L' and <r|sL> */
  dfdist_.resize(L.size());
  for (int i = 0; i != L.size(); ++i) {
    vector<shared_ptr<const Atom>> atoms1;
    atoms1.resize(atoms0.size());
    int iat = 0;
    for (auto& atom : atoms0) {
      atoms1[i] = make_shared<const Atom>(*atom, lattice_vectors_[i]);
      ++iat;
    }
    dfdist_[i] = make_shared<PDFDist_ints>(L, nbas, naux, atoms0, atoms1, aux_atoms, thresh, data1_);
  }
}


void PDFDist::compute_aux_charge(const vector<shared_ptr<const Shell>>& ashell) {

  TaskQueue<PDFIntTask_aux> tasks(ashell.size());
  data1_ = make_shared<VectorB>(naux_);

  // <a|.>
  auto i1 = make_shared<const Shell>(ashell.front()->spherical());
  int j0 = 0;
  for (auto& i0 : ashell) {
    tasks.emplace_back(array<shared_ptr<const Shell>,2>{{i1, i0}}, j0, data1_);
    j0 += i0->nbasis();
  }

  tasks.compute();

  const double q = data1_->rms() * naux_;
  for (auto& idata : *data1_) idata /= q;

  // P_{ij} = <i|.><.|j>
  Matrix p(naux_, naux_);
  dger_(naux_, naux_, 1.0, data1_->data(), 1, data1_->data(), 1, p.data(), naux_);
  projector_ = make_shared<const Matrix>(p);
}


void PDFDist::pcompute_2index(const vector<shared_ptr<const Shell>>& ashell, const double throverlap) {

  Timer time;

  TaskQueue<PDFIntTask_2index> tasks(ashell.size() * ashell.size() * ncell());

  data2_ = make_shared<Matrix>(naux_, naux_, serial_);
  auto b3 = make_shared<const Shell>(ashell.front()->spherical());

  int u = 0;
  int o0 = 0;
  for (auto& b0 : ashell) {
    int o1 = 0;
    for (auto& b1 : ashell) {
      for (auto& L : lattice_vectors_) {
        auto b11 = make_shared<const Shell>(*(b1->move_atom(L)));
        if ((u++ % mpi__->size() == mpi__->rank()) || serial_)
          tasks.emplace_back(array<shared_ptr<const Shell>,4>{{b11, b3, b0, b3}}, array<int,2>{{o0, o1}}, data2_);
      }
      o1 += b1->nbasis();
    }
    o0 += b0->nbasis();
  }

  time.tick_print("2-index integrals prep");
  tasks.compute();

  if (!serial_)
    data2_->allreduce();

  time.tick_print("2-index integrals");

  if (!projector_)
    throw logic_error("failed attempt to project data2_ before computing the projection matrix");
  *data2_ = *projector_ * *data2_ * *projector_;

  // P_C = 1 - P
  auto projectorC = make_shared<Matrix>(*projector_);
  projectorC->unit();
  *projectorC = *projectorC - *projector_;

  // make data2_ positive definite
  *data2_ += *projectorC;
//data2_->inverse_half(throverlap);
  data2_->inverse(); //TODO: linear dependency can be a problem!
  // use data2_ within node
  data2_->localize();
  time.tick_print("computing inverse");
}
