//
// BAGEL - Parallel electron correlation program.
// Filename: gradtask.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/grad/gradeval_base.h>
#include <src/integral/rys/gradbatch.h>
#include <src/integral/rys/gnaibatch.h>
#include <src/integral/rys/gsmallnaibatch.h>
#include <src/integral/rys/gsmalleribatch.h>
#include <src/integral/os/goverlapbatch.h>
#include <src/integral/os/gkineticbatch.h>
#include <src/smith/prim_op.h>
#ifdef LIBINT_INTERFACE
  #include <src/integral/libint/glibint.h>
#endif

using namespace std;
using namespace bagel;


void GradTask3::compute() {
#ifdef LIBINT_INTERFACE
  GLibint gradbatch(shell_);
#else
  GradBatch gradbatch(shell_, 0.0);
#endif
  gradbatch.compute();
  const size_t block = gradbatch.size_block();
  const size_t sblock = shell_[1]->nbasis()*shell_[2]->nbasis()*shell_[3]->nbasis();
  assert(sblock <= block);

  // unfortunately the convention is different...
  array<int,4> jatom = {{-1, atomindex_[2], atomindex_[1], atomindex_[0]}};
  if (gradbatch.swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
  if (gradbatch.swap01()) swap(jatom[0], jatom[1]);
  if (gradbatch.swap23()) swap(jatom[2], jatom[3]);

  shared_ptr<Matrix> db1 = den_->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis());
  shared_ptr<Matrix> db2 = den_->get_block(offset_[2], shell_[1]->nbasis(), offset_[0], shell_[3]->nbasis(), offset_[1], shell_[2]->nbasis());
  SMITH::sort_indices<0,2,1,1,1,1,1>(db2->data(), db1->data(), shell_[1]->nbasis(), shell_[3]->nbasis(), shell_[2]->nbasis());

  for (int iatom = 0; iatom != 4; ++iatom) {
    if (jatom[iatom] < 0) continue;
    array<double,3> sum = {{0.0, 0.0, 0.0}};
    for (int icart = 0; icart != 3; ++icart) {
      const double* ppt = gradbatch.data(icart+iatom*3);
      sum[icart] += inner_product(ppt, ppt+sblock, db1->data(), 0.0);
    }
    lock_guard<mutex> lock(ge_->mutex_[jatom[iatom]]);
    for (int icart = 0; icart != 3; ++icart)
      ge_->grad_->element(icart, jatom[iatom]) += 0.5 * sum[icart] * (shell_[2] == shell_[3] ? 1.0 : 2.0);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// significant repetition, but not sure how to nicely eliminate them
void GradTask1f::compute() {
#ifdef LIBINT_INTERFACE
  GLibint gradbatch(shell_);
#else
  GradBatch gradbatch(shell_, 0.0);
#endif
  gradbatch.compute();
  const size_t block = gradbatch.size_block();
  const size_t sblock = shell_[2]->nbasis()*shell_[3]->nbasis();
  assert(sblock <= block && shell_[1]->nbasis() == 1);

  // unfortunately the convention is different...
  array<int,4> jatom = {{-1, atomindex_[2], atomindex_[1], atomindex_[0]}};
  if (gradbatch.swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
  if (gradbatch.swap01()) swap(jatom[0], jatom[1]);
  if (gradbatch.swap23()) swap(jatom[2], jatom[3]);

  shared_ptr<Matrix> db1 = den_->get_submatrix(offset_[1], offset_[0], shell_[2]->nbasis(), shell_[3]->nbasis());

  for (int iatom = 0; iatom != 4; ++iatom) {
    if (jatom[iatom] < 0) continue;
    array<double,3> sum = {{0.0, 0.0, 0.0}};
    for (int icart = 0; icart != 3; ++icart) {
      const double* ppt = gradbatch.data(icart+iatom*3);
      sum[icart] += inner_product(ppt, ppt+sblock, db1->data(), 0.0);
    }
    lock_guard<mutex> lock(ge_->mutex_[jatom[iatom]]);
    for (int icart = 0; icart != 3; ++icart)
      ge_->grad_->element(icart, jatom[iatom]) += sum[icart];
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GradTask2::compute() {
#ifdef LIBINT_INTERFACE
  GLibint gradbatch(shell_);
#else
  GradBatch gradbatch(shell_, 0.0);
#endif
  gradbatch.compute();
  const size_t block = gradbatch.size_block();

  // unfortunately the convention is different...
  int jatom[4] = {atomindex_[1], -1, atomindex_[0], -1};
  if (gradbatch.swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
  if (gradbatch.swap01()) swap(jatom[0], jatom[1]);
  if (gradbatch.swap23()) swap(jatom[2], jatom[3]);

  for (int iatom = 0; iatom != 4; ++iatom) {
    if (jatom[iatom] < 0) continue;
    array<double,3> sum = {{0.0, 0.0, 0.0}};
    for (int icart = 0; icart != 3; ++icart) {
      const double* ppt = gradbatch.data(icart+iatom*3);
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[2]->nbasis(); ++j0) {
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++ppt) {
          sum[icart] += *ppt * den2_->element(j1,j0);
          sum[icart] += *ppt * den2_->element(j0,j1);
        }
      }
    }
    lock_guard<mutex> lock(ge_->mutex_[jatom[iatom]]);
    // first 0.5 from symmetrization. second 0.5 from the Hamiltonian
    for (int icart = 0; icart != 3; ++icart)
      ge_->grad_->element(icart, jatom[iatom]) -= 0.5 * sum[icart] * 0.5 * (shell_[0] == shell_[2] ? 1.0 : 2.0);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GradTask1::compute() {
  auto grad_local = make_shared<GradFile>(ge_->geom_->natom());
  *grad_local += *compute_nai();
  *grad_local += *compute_os<GKineticBatch>(den3_);
  *grad_local -= *compute_os<GOverlapBatch>(eden_);

  for (int iatom = 0; iatom != ge_->geom_->natom(); ++iatom) {
    lock_guard<mutex> lock(ge_->mutex_[iatom]);
    ge_->grad_->element(0, iatom) += grad_local->element(0, iatom);
    ge_->grad_->element(1, iatom) += grad_local->element(1, iatom);
    ge_->grad_->element(2, iatom) += grad_local->element(2, iatom);
  }
}


shared_ptr<GradFile> GradTask1::compute_nai() const {
  const int dimb1 = shell_[0]->nbasis();
  const int dimb0 = shell_[1]->nbasis();
  GNAIBatch batch2(shell_, ge_->geom_, tie(atomindex_[1], atomindex_[0]));
  batch2.compute();
  shared_ptr<Matrix> cden = den2_->get_submatrix(offset_[1], offset_[0], dimb1, dimb0);
  const int dummy = -1;
  return batch2.compute_gradient(cden, dummy, dummy, ge_->geom_->natom());
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GradTask3r::compute() {
  shared_ptr<GradFile> grad_local = compute_smalleri();
  list<int> done;
  for (int i = 0; i != 3; ++i) {
    const int iatom = atomindex_[i];
    if (find(done.begin(), done.end(), iatom) != done.end()) continue; // should not add twice
    done.push_back(iatom);

    lock_guard<mutex> lock(ge_->mutex_[iatom]);
    ge_->grad_->element(0, iatom) += grad_local->element(0, iatom);
    ge_->grad_->element(1, iatom) += grad_local->element(1, iatom);
    ge_->grad_->element(2, iatom) += grad_local->element(2, iatom);
  }
}


shared_ptr<GradFile> GradTask3r::compute_smalleri() const {
  GSmallERIBatch batch(shell_, array<int,3>{{atomindex_[0], atomindex_[1], atomindex_[2]}}, ge_->geom_->natom());
  batch.compute();
  array<shared_ptr<const Matrix>,6> d = {{
    rden3_[0]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()),
    rden3_[1]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()),
    rden3_[2]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()),
    rden3_[3]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()),
    rden3_[4]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()),
    rden3_[5]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()) }};
  return batch.compute_gradient(d);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GradTask1rf::compute() {
  shared_ptr<GradFile> grad_local = compute_smalleri();
  list<int> done;
  for (int i = 0; i != 3; ++i) {
    const int iatom = atomindex_[i];
    if (find(done.begin(), done.end(), iatom) != done.end()) continue; // should not add twice
    done.push_back(iatom);

    lock_guard<mutex> lock(ge_->mutex_[iatom]);
    ge_->grad_->element(0, iatom) += grad_local->element(0, iatom);
    ge_->grad_->element(1, iatom) += grad_local->element(1, iatom);
    ge_->grad_->element(2, iatom) += grad_local->element(2, iatom);
  }
}


shared_ptr<GradFile> GradTask1rf::compute_smalleri() const {
  GSmallERIBatch batch(shell_, array<int,3>{{atomindex_[0], atomindex_[1], atomindex_[2]}}, ge_->geom_->natom());
  batch.compute();
  array<shared_ptr<const Matrix>,6> d = {{
    rden_[0]->get_submatrix(offset_[1], offset_[0], shell_[2]->nbasis(), shell_[3]->nbasis()),
    rden_[1]->get_submatrix(offset_[1], offset_[0], shell_[2]->nbasis(), shell_[3]->nbasis()),
    rden_[2]->get_submatrix(offset_[1], offset_[0], shell_[2]->nbasis(), shell_[3]->nbasis()),
    rden_[3]->get_submatrix(offset_[1], offset_[0], shell_[2]->nbasis(), shell_[3]->nbasis()),
    rden_[4]->get_submatrix(offset_[1], offset_[0], shell_[2]->nbasis(), shell_[3]->nbasis()),
    rden_[5]->get_submatrix(offset_[1], offset_[0], shell_[2]->nbasis(), shell_[3]->nbasis()) }};
  return batch.compute_gradient(d);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GradTask1r::compute() {
  shared_ptr<GradFile> grad_local = compute_smallnai();
  for (int iatom = 0; iatom != ge_->geom_->natom(); ++iatom) {
    lock_guard<mutex> lock(ge_->mutex_[iatom]);
    ge_->grad_->element(0, iatom) += grad_local->element(0, iatom);
    ge_->grad_->element(1, iatom) += grad_local->element(1, iatom);
    ge_->grad_->element(2, iatom) += grad_local->element(2, iatom);
  }
}


shared_ptr<GradFile> GradTask1r::compute_smallnai() const {
  const int dimb1 = shell_[0]->nbasis();
  const int dimb0 = shell_[1]->nbasis();
  GSmallNAIBatch batch(shell_, ge_->geom_, tie(atomindex_[1], atomindex_[0]));
  batch.compute();

  array<shared_ptr<const Matrix>,6> dmat;
  auto iter = rden_.begin();
  for (auto& i : dmat) {
    shared_ptr<Matrix> tmp = (*iter)->get_submatrix(offset_[1], offset_[0], dimb1, dimb0);
    tmp->localize();
    i = tmp;
    ++iter;
  }
  return batch.compute_gradient(dmat);
}
