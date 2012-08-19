//
// Newint - Parallel electron correlation program.
// Filename: gradtask.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <array>
#include <src/grad/gradeval_base.h>

using namespace std;

void GradTask::compute() {
  if (atomindex_.size() == 3) {
    // pointer to stack
    GradBatch gradbatch(shell_, 0.0);
    gradbatch.compute();
    const size_t block = gradbatch.size_block();

    // unfortunately the convention is different...
    array<int,4> jatom = {{-1, atomindex_[2], atomindex_[1], atomindex_[0]}};
    if (gradbatch.swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
    if (gradbatch.swap01()) swap(jatom[0], jatom[1]);
    if (gradbatch.swap23()) swap(jatom[2], jatom[3]);

    for (int i = 0; i != 12; ++i) {
      // if this is a dummy atom
      if (jatom[i/3] < 0) continue;

      const double* ppt = gradbatch.data() + i*block;
      double sum = 0.0;
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0) {  
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1) {  
          for (int j2 = offset_[2]; j2 != offset_[2] + shell_[1]->nbasis(); ++j2, ++ppt) {  
            // first we need to have a scheme to receive blocks before accessing the elements
            sum += *ppt * *den_->ptr(j2, j1, j0);
            sum += *ppt * *den_->ptr(j2, j0, j1);
          }
        }
      }
      boost::lock_guard<boost::mutex> lock(ge_->mutex_[jatom[i/3]]);
      ge_->grad_->data(jatom[i/3], i%3) += 0.5 * sum * (shell_[2] == shell_[3] ? 1.0 : 2.0);
    }
  } else if (atomindex_.size() == 2) {
    // pointer to stack
    GradBatch gradbatch(shell_, 0.0);
    gradbatch.compute();
    const size_t block = gradbatch.size_block();

    // unfortunately the convention is different...
    int jatom[4] = {atomindex_[1], -1, atomindex_[0], -1};
    if (gradbatch.swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
    if (gradbatch.swap01()) swap(jatom[0], jatom[1]);
    if (gradbatch.swap23()) swap(jatom[2], jatom[3]);

    for (int i = 0; i != 12; ++i) {
      // if this is a dummy atom
      if (jatom[i/3] < 0) continue;

      const double* ppt = gradbatch.data() + i*block;
      double sum = 0.0;
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[2]->nbasis(); ++j0) {  
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++ppt) {  
          sum += *ppt * den2_->element(j1,j0);
          sum += *ppt * den2_->element(j0,j1);
        }
      }
      // first 0.5 from symmetrization. second 0.5 from the Hamiltonian
      boost::lock_guard<boost::mutex> lock(ge_->mutex_[jatom[i/3]]);
      ge_->grad_->data(jatom[i/3],i%3) -= 0.5 * sum * 0.5 * (shell_[0] == shell_[2] ? 1.0 : 2.0);
    }
  } else {
    throw logic_error("calling GradTask::compute() with illegal setups");
  }
} 

