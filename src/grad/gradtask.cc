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


#include <cstddef>
#include <array>
#include <src/grad/gradeval_base.h>
#include <src/grad/gradfile.h>
#include <src/grad/gnaibatch.h>
#include <src/grad/goverlapbatch.h>
#include <src/grad/gkineticbatch.h>
#ifdef LIBINT_INTERFACE
  #include <src/grad/glibint.h>
#endif

using namespace std;
using namespace bagel;

void GradTask::compute() {
  if (rank_ == 1) {
    const int iatom0 = atomindex_[0];
    const int iatom1 = atomindex_[1];
    const int nbasis = ge_->geom_->nbasis();
    shared_ptr<GradFile> grad_local(new GradFile(ge_->geom_->natom()));
    const int dimb1 = shell2_[0]->nbasis();
    const int dimb0 = shell2_[1]->nbasis();
    {
      GNAIBatch batch2(shell2_, ge_->geom_, tie(iatom1, iatom0));
      batch2.compute();
      const double* ndata = batch2.data();
      const size_t s = batch2.size_block();
      for (int ia = 0; ia != ge_->geom_->natom()*3; ++ia) {
        for (int i = offset_[0], cnt = 0; i != dimb0 + offset_[0]; ++i) {
          for (int j = offset_[1]; j != dimb1 + offset_[1]; ++j, ++cnt) {
            grad_local->data(ia) += ndata[cnt+s*ia] * den2_->data(i*nbasis+j);
          }
        }
      }
    }
    {
      GKineticBatch batch(shell2_, ge_->geom_);
      const double* kdata = batch.data();
      batch.compute();
      const size_t s = batch.size_block();
      for (int i = offset_[0], cnt = 0; i != dimb0 + offset_[0]; ++i) {
        for (int j = offset_[1]; j != dimb1 + offset_[1]; ++j, ++cnt) {
          int jatom0 = batch.swap01() ? iatom1 : iatom0;
          int jatom1 = batch.swap01() ? iatom0 : iatom1;
          for (int k = 0; k != 3; ++k) {
            grad_local->data(3*jatom1+k) += kdata[cnt+s*k    ] * den2_->data(i*nbasis+j);
            grad_local->data(3*jatom0+k) += kdata[cnt+s*(k+3)] * den2_->data(i*nbasis+j);
          }
        }
      }
    }
    {
      GOverlapBatch batch(shell2_, ge_->geom_);
      const double* odata = batch.data();
      batch.compute();
      const size_t s = batch.size_block();
      for (int i = offset_[0], cnt = 0; i != dimb0 + offset_[0]; ++i) {
        for (int j = offset_[1]; j != dimb1 + offset_[1]; ++j, ++cnt) {
          int jatom0 = batch.swap01() ? iatom1 : iatom0;
          int jatom1 = batch.swap01() ? iatom0 : iatom1;
          for (int k = 0; k != 3; ++k) {
            grad_local->data(3*jatom1+k) -= odata[cnt+s*k    ] * eden_->data(i*nbasis+j);
            grad_local->data(3*jatom0+k) -= odata[cnt+s*(k+3)] * eden_->data(i*nbasis+j);
          }
        }
      }
    }
    for (int iatom = 0; iatom != ge_->geom_->natom(); ++iatom) {
      boost::lock_guard<boost::mutex> lock(ge_->mutex_[iatom]);
      ge_->grad_->data(3*iatom+0) += grad_local->data(3*iatom+0);
      ge_->grad_->data(3*iatom+1) += grad_local->data(3*iatom+1);
      ge_->grad_->data(3*iatom+2) += grad_local->data(3*iatom+2);
    }

  } else if (rank_ == 3) {

#ifdef LIBINT_INTERFACE
    GLibint gradbatch(shell_);
#else
    GradBatch gradbatch(shell_, 0.0);
#endif
    gradbatch.compute();
    const size_t block = gradbatch.size_block();


    // unfortunately the convention is different...
    array<int,4> jatom = {{-1, atomindex_[2], atomindex_[1], atomindex_[0]}};
    if (gradbatch.swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
    if (gradbatch.swap01()) swap(jatom[0], jatom[1]);
    if (gradbatch.swap23()) swap(jatom[2], jatom[3]);

    for (int iatom = 0; iatom != 4; ++iatom) {
      if (jatom[iatom] < 0) continue;
      array<double,3> sum = {{0.0, 0.0, 0.0}};
      for (int icart = 0; icart != 3; ++icart) {
        const double* ppt = gradbatch.data() + (icart+iatom*3)*block;
        for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0) {
          for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1) {
            for (int j2 = offset_[2]; j2 != offset_[2] + shell_[1]->nbasis(); ++j2, ++ppt) {
              // first we need to have a scheme to receive blocks before accessing the elements
              sum[icart] += *ppt * *den_->ptr(j2, j1, j0);
              sum[icart] += *ppt * *den_->ptr(j2, j0, j1);
            }
          }
        }
      }
      boost::lock_guard<boost::mutex> lock(ge_->mutex_[jatom[iatom]]);
      for (int icart = 0; icart != 3; ++icart)
        ge_->grad_->data(jatom[iatom], icart) += 0.5 * sum[icart] * (shell_[2] == shell_[3] ? 1.0 : 2.0);
    }
  } else if (rank_ == 2) {
    // pointer to stack
    GradBatch gradbatch(shell_, 0.0);
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
      boost::lock_guard<boost::mutex> lock(ge_->mutex_[jatom[iatom]]);
      // first 0.5 from symmetrization. second 0.5 from the Hamiltonian
      for (int icart = 0; icart != 3; ++icart)
        ge_->grad_->data(jatom[iatom],icart) -= 0.5 * sum[icart] * 0.5 * (shell_[0] == shell_[2] ? 1.0 : 2.0);
    }
  } else {
    throw logic_error("calling GradTask::compute() with illegal setups");
  }
}

