//
// Newint - Parallel electron correlation program.
// Filename: gradeval_base.cc
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


#include <src/grad/gradeval_base.h>
#include <src/grad/gnaibatch.h>
#include <src/grad/goverlapbatch.h>
#include <src/grad/gkineticbatch.h>

using namespace std;

void GradEval_base::compute_grad1e() {

  const int natom = geom_->natom();

  const vector<shared_ptr<Atom> > atoms = geom_->atoms(); 
  const vector<vector<int> > offsets = geom_->offsets();
  const int nbasis = geom_->nbasis();

  // only lower half will be stored
  for (int iatom0 = 0; iatom0 != natom; ++iatom0) {
    // iatom1 = iatom1;
    const shared_ptr<Atom> catom0 = atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = offsets[iatom0];
    const vector<shared_ptr<Shell> > shell0 = catom0->shells();

    for (int iatom1 = 0; iatom1 != natom; ++iatom1) {
      const shared_ptr<Atom> catom1 = atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = offsets[iatom1];
      const vector<shared_ptr<Shell> > shell1 = catom1->shells();

      for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
        const int offset0 = coffset0[ibatch0]; 
        shared_ptr<Shell> b0 = shell0[ibatch0];
        for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
          const int offset1 = coffset1[ibatch1]; 
          shared_ptr<Shell> b1 = shell1[ibatch1];
          vector<shared_ptr<Shell> > input;
          input.push_back(b1);
          input.push_back(b0);

          const int dimb1 = input[0]->nbasis(); 
          const int dimb0 = input[1]->nbasis(); 
          const int ang1 = b1->angular_number();
          const int ang0 = b0->angular_number();
          const int s = (ang1+1)*(ang1+2)*(ang0+1)*(ang0+2)/4 * b1->num_primitive()*b0->num_primitive();

          {
            GNAIBatch batch2(input, geom_, tie(iatom1, iatom0));
            batch2.compute();
            const double* ndata = batch2.data();
            for (int ia = 0; ia != natom*3; ++ia) {
              for (int i = offset0, cnt = 0; i != dimb0 + offset0; ++i) {
                for (int j = offset1; j != dimb1 + offset1; ++j, ++cnt) {
                  grad1e(ia)->data(i*nbasis+j) += ndata[cnt+s*(ia)];
                }
              }
            }
          }
          {
            GKineticBatch batch(input, geom_);
            const double* kdata = batch.data();
            batch.compute();
            for (int i = offset0, cnt = 0; i != dimb0 + offset0; ++i) {
              for (int j = offset1; j != dimb1 + offset1; ++j, ++cnt) {
                int jatom0 = batch.swap01() ? iatom1 : iatom0;
                int jatom1 = batch.swap01() ? iatom0 : iatom1;
                grad1e(3*jatom1+0)->data(i*nbasis+j) += kdata[cnt];
                grad1e(3*jatom1+1)->data(i*nbasis+j) += kdata[cnt+s];
                grad1e(3*jatom1+2)->data(i*nbasis+j) += kdata[cnt+s*2];
                grad1e(3*jatom0+0)->data(i*nbasis+j) += kdata[cnt+s*3];
                grad1e(3*jatom0+1)->data(i*nbasis+j) += kdata[cnt+s*4];
                grad1e(3*jatom0+2)->data(i*nbasis+j) += kdata[cnt+s*5];
              }
            }
          }
          {
            GOverlapBatch batch(input, geom_);
            const double* odata = batch.data();
            batch.compute();
            for (int i = offset0, cnt = 0; i != dimb0 + offset0; ++i) {
              for (int j = offset1; j != dimb1 + offset1; ++j, ++cnt) {
                int jatom0 = batch.swap01() ? iatom1 : iatom0;
                int jatom1 = batch.swap01() ? iatom0 : iatom1;
                grado(3*jatom1+0)->data(i*nbasis+j) += odata[cnt];
                grado(3*jatom1+1)->data(i*nbasis+j) += odata[cnt+s];
                grado(3*jatom1+2)->data(i*nbasis+j) += odata[cnt+s*2];
                grado(3*jatom0+0)->data(i*nbasis+j) += odata[cnt+s*3];
                grado(3*jatom0+1)->data(i*nbasis+j) += odata[cnt+s*4];
                grado(3*jatom0+2)->data(i*nbasis+j) += odata[cnt+s*5];
              }
            }
          }

        }
      } 
    }
  }
}
