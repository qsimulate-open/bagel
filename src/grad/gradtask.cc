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


#include <stddef.h>
#include <array>
#include <src/grad/gradeval_base.h>
#include <src/grad/gradfile.h>
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


shared_ptr<GradFile> GradTask::compute_nai() const {
  const int dimb1 = shell2_[0]->nbasis();
  const int dimb0 = shell2_[1]->nbasis();
  GNAIBatch batch2(shell2_, ge_->geom_, tie(atomindex_[1], atomindex_[0]));
  batch2.compute();
  shared_ptr<Matrix> cden = den2_->get_submatrix(offset_[1], offset_[0], dimb1, dimb0);
  const int dummy = -1;
  return batch2.compute_gradient(cden, dummy, dummy, ge_->geom_->natom());
}


shared_ptr<GradFile> GradTask::compute_smallnai() const {
  const int dimb1 = shell2_[0]->nbasis();
  const int dimb0 = shell2_[1]->nbasis();
  GSmallNAIBatch batch(shell2_, ge_->geom_, tie(atomindex_[1], atomindex_[0]));
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


shared_ptr<GradFile> GradTask::compute_smalleri() const {
  GSmallERIBatch batch(shell_, array<int,3>{{atomindex_[0], atomindex_[1], atomindex_[2]}}, ge_->geom_->natom());
  batch.compute();

  array<unique_ptr<double[]>,6> d = {{
    rden3_[0]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()),
    rden3_[1]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()),
    rden3_[2]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()),
    rden3_[3]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()),
    rden3_[4]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()),
    rden3_[5]->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis()) }};

  return batch.compute_gradient(d);
}


void GradTask::compute() {

  if (rank_ == 1) {
    auto grad_local = make_shared<GradFile>(ge_->geom_->natom());
    *grad_local += *compute_nai();
    *grad_local += *compute_os<GKineticBatch>(den3_);
    *grad_local -= *compute_os<GOverlapBatch>(eden_);

    for (int iatom = 0; iatom != ge_->geom_->natom(); ++iatom) {
      lock_guard<mutex> lock(ge_->mutex_[iatom]);
      ge_->grad_->data(0, iatom) += grad_local->data(0, iatom);
      ge_->grad_->data(1, iatom) += grad_local->data(1, iatom);
      ge_->grad_->data(2, iatom) += grad_local->data(2, iatom);
    }

  // relativistic one-electron integrals
  } else if (rank_ == -1) {

    shared_ptr<GradFile> grad_local = compute_smallnai();
    for (int iatom = 0; iatom != ge_->geom_->natom(); ++iatom) {
      lock_guard<mutex> lock(ge_->mutex_[iatom]);
      ge_->grad_->data(0, iatom) += grad_local->data(0, iatom);
      ge_->grad_->data(1, iatom) += grad_local->data(1, iatom);
      ge_->grad_->data(2, iatom) += grad_local->data(2, iatom);
    }

  } else if (rank_ == -3) {

    shared_ptr<GradFile> grad_local = compute_smalleri();
    list<int> done;
    for (int i = 0; i != 3; ++i) {
      const int iatom = atomindex_[i];
      if (find(done.begin(), done.end(), iatom) != done.end()) continue; // should not add twice
      done.push_back(iatom);

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

    unique_ptr<double[]> db1 = den_->get_block(offset_[2], shell_[1]->nbasis(), offset_[1], shell_[2]->nbasis(), offset_[0], shell_[3]->nbasis());
    unique_ptr<double[]> db2 = den_->get_block(offset_[2], shell_[1]->nbasis(), offset_[0], shell_[3]->nbasis(), offset_[1], shell_[2]->nbasis());
    SMITH::sort_indices<0,2,1,1,1,1,1>(db2, db1, shell_[1]->nbasis(), shell_[3]->nbasis(), shell_[2]->nbasis());

    for (int iatom = 0; iatom != 4; ++iatom) {
      if (jatom[iatom] < 0) continue;
      array<double,3> sum = {{0.0, 0.0, 0.0}};
      for (int icart = 0; icart != 3; ++icart) {
        const double* ppt = gradbatch.data(icart+iatom*3);
        sum[icart] += ddot_(sblock, ppt, 1, db1.get(), 1);
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
        ge_->grad_->data(icart, jatom[iatom]) -= 0.5 * sum[icart] * 0.5 * (shell_[0] == shell_[2] ? 1.0 : 2.0);
    }
  } else {
    throw logic_error("calling GradTask::compute() with illegal setups");
  }
}

