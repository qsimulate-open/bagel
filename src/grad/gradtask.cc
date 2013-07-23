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
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <stddef.h>
#include <array>
#include <src/grad/gradeval_base.h>
#include <src/grad/gradfile.h>
#include <src/integral/rys/gradbatch.h>
#include <src/integral/rys/gnaibatch.h>
#include <src/integral/os/goverlapbatch.h>
#include <src/integral/os/gkineticbatch.h>
#include <src/smith/prim_op.h>
#ifdef LIBINT_INTERFACE
  #include <src/integral/libint/glibint.h>
#endif

using namespace std;
using namespace bagel;


shared_ptr<GradFile> GradTask::compute_nai(shared_ptr<const Matrix> den) const {
  const int iatom0 = atomindex_[0];
  const int iatom1 = atomindex_[1];
  const int nbasis = ge_->geom_->nbasis();
  auto grad_local = make_shared<GradFile>(ge_->geom_->natom());
  const int dimb1 = shell2_[0]->nbasis();
  const int dimb0 = shell2_[1]->nbasis();

  GNAIBatch batch2(shell2_, ge_->geom_, tie(iatom1, iatom0));
  batch2.compute();

  const double* ndata = batch2.data();
  const size_t s = batch2.size_block();
  for (int ia = 0; ia != ge_->geom_->natom()*3; ++ia) {
    for (int i = offset_[0], cnt = 0; i != dimb0 + offset_[0]; ++i) {
      for (int j = offset_[1]; j != dimb1 + offset_[1]; ++j, ++cnt) {
        grad_local->data(ia) += ndata[cnt+s*ia] * den->data(i*nbasis+j);
      }
    }
  }
  return grad_local;
}


void GradTask::compute() {

  if (rank_ == 1) {
    auto grad_local = make_shared<GradFile>(ge_->geom_->natom());
    *grad_local += *compute_nai(den2_);
    *grad_local += *compute_os<GKineticBatch>(den3_);
    *grad_local -= *compute_os<GOverlapBatch>(eden_);

    for (int iatom = 0; iatom != ge_->geom_->natom(); ++iatom) {
      lock_guard<mutex> lock(ge_->mutex_[iatom]);
      ge_->grad_->data(0, iatom) += grad_local->data(0, iatom);
      ge_->grad_->data(1, iatom) += grad_local->data(1, iatom);
      ge_->grad_->data(2, iatom) += grad_local->data(2, iatom);
    }

  } else if (rank_ == 3) {
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

    const unique_ptr<double[]> db1 = den_->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis());
    const unique_ptr<double[]> db2 = den_->get_block(offset_[2], shell_[1]->nbasis(), offset_[0], shell_[3]->nbasis(), offset_[1], shell_[2]->nbasis());
    unique_ptr<double[]> db3(new double[sblock]);
    SMITH::sort_indices<0,2,1,0,1,1,1>(db2, db3, shell_[1]->nbasis(), shell_[3]->nbasis(), shell_[2]->nbasis());

    for (int iatom = 0; iatom != 4; ++iatom) {
      if (jatom[iatom] < 0) continue;
      array<double,3> sum = {{0.0, 0.0, 0.0}};
      for (int icart = 0; icart != 3; ++icart) {
        const double* ppt = gradbatch.data() + (icart+iatom*3)*block;
        sum[icart] += ddot_(sblock, ppt, 1, db1.get(), 1);
        sum[icart] += ddot_(sblock, ppt, 1, db3.get(), 1);
      }
      lock_guard<mutex> lock(ge_->mutex_[jatom[iatom]]);
      for (int icart = 0; icart != 3; ++icart)
        ge_->grad_->data(icart, jatom[iatom]) += 0.5 * sum[icart] * (shell_[2] == shell_[3] ? 1.0 : 2.0);
    }
  } else if (rank_ == 2) {
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
        const double* ppt = gradbatch.data() + (icart+iatom*3)*block;
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
        ge_->grad_->data(icart, jatom[iatom]) -= 0.5 * sum[icart] * 0.5 * (shell_[0] == shell_[2] ? 1.0 : 2.0);
    }
  } else {
    throw logic_error("calling GradTask::compute() with illegal setups");
  }
}

