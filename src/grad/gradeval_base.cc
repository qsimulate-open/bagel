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


shared_ptr<GradFile> GradEval_base::contract_gradient(const shared_ptr<const Matrix1e> d, const shared_ptr<const Matrix1e> w,
                                                      const shared_ptr<const DF_AO> o, const unique_ptr<double[]>& o2) const {

  shared_ptr<Grad1eFile> g1(new Grad1eFile(geom_));
  shared_ptr<Grad1eFile> go(new Grad1eFile(geom_));
  compute_grad1e_integrals(g1, go);

  vector<double> grad = contract_grad1e(d, w, g1, go); 
  vector<double> tmp0 = contract_grad2e(o);
  vector<double> tmp1 = contract_grad2e_2index(o2);

  for (auto i = grad.begin(), t0 = tmp0.begin(), t1 = tmp1.begin(); i != grad.end(); ++i, ++t0, ++t1)
    *i += *t0 + *t1;

  return shared_ptr<GradFile>(new GradFile(grad));
}


void GradEval_base::compute_grad1e_integrals(shared_ptr<Grad1eFile> g1, shared_ptr<Grad1eFile> go) const {

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
                  g1->data(ia)->data(i*nbasis+j) += ndata[cnt+s*(ia)];
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
                g1->data(3*jatom1+0)->data(i*nbasis+j) += kdata[cnt];
                g1->data(3*jatom1+1)->data(i*nbasis+j) += kdata[cnt+s];
                g1->data(3*jatom1+2)->data(i*nbasis+j) += kdata[cnt+s*2];
                g1->data(3*jatom0+0)->data(i*nbasis+j) += kdata[cnt+s*3];
                g1->data(3*jatom0+1)->data(i*nbasis+j) += kdata[cnt+s*4];
                g1->data(3*jatom0+2)->data(i*nbasis+j) += kdata[cnt+s*5];
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
                go->data(3*jatom1+0)->data(i*nbasis+j) += odata[cnt];
                go->data(3*jatom1+1)->data(i*nbasis+j) += odata[cnt+s];
                go->data(3*jatom1+2)->data(i*nbasis+j) += odata[cnt+s*2];
                go->data(3*jatom0+0)->data(i*nbasis+j) += odata[cnt+s*3];
                go->data(3*jatom0+1)->data(i*nbasis+j) += odata[cnt+s*4];
                go->data(3*jatom0+2)->data(i*nbasis+j) += odata[cnt+s*5];
              }
            }
          }

        }
      } 
    }
  }
}

vector<double> GradEval_base::contract_grad1e(const shared_ptr<const Matrix1e> d, const shared_ptr<const Matrix1e> w,
                                              const shared_ptr<const Grad1eFile> g1, const shared_ptr<const Grad1eFile> go) const {

  // first Vnuc derivative
  vector<double> grad = geom_->compute_grad_vnuc();
  assert(grad.size() == geom_->natom()*3);

  int i = 0;
  for (auto g = grad.begin(); g != grad.end(); ++g, ++i)
    *g += g1->data(i)->ddot(d) - go->data(i)->ddot(w);

  return grad;
}


vector<double> GradEval_base::contract_grad2e(const shared_ptr<const DF_AO> o) const {

  vector<double> grad(geom_->natom()*3, 0.0);

  // loop over atoms (using symmetry b0 <-> b1)
  vector<shared_ptr<Atom> > aux_atoms = geom_->aux_atoms();
  vector<vector<int> > aux_offsets = geom_->aux_offsets();
  const vector<shared_ptr<Atom> > atoms = geom_->atoms(); 
  const vector<vector<int> > offsets = geom_->offsets();
  const int nbasis = geom_->nbasis();
  for (int iatom0 = 0; iatom0 != geom_->natom(); ++iatom0) {
    const shared_ptr<Atom> catom0 = atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = offsets[iatom0];
    const vector<shared_ptr<Shell> > shell0 = catom0->shells();

    for (int iatom1 = iatom0; iatom1 != geom_->natom(); ++iatom1) {
      const shared_ptr<Atom> catom1 = atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = offsets[iatom1];
      const vector<shared_ptr<Shell> > shell1 = catom1->shells();

      for (int iatom2 = 0; iatom2 != geom_->natom(); ++iatom2) {
        const shared_ptr<Atom> catom2 = aux_atoms[iatom2];
        const int numshell2 = catom2->shells().size();
        const vector<int> coffset2 = aux_offsets[iatom2];
        const vector<shared_ptr<Shell> > shell2 = catom2->shells();


        // dummy shell
        const shared_ptr<Shell> b3(new Shell(shell2.front()->spherical()));

        for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
          const int offset0 = coffset0[ibatch0]; 
          shared_ptr<Shell> b0 = shell0[ibatch0];

          for (int ibatch1 = (iatom0 == iatom1 ? ibatch0 : 0); ibatch1 != numshell1; ++ibatch1) {
            const int offset1 = coffset1[ibatch1]; 
            shared_ptr<Shell> b1 = shell1[ibatch1];

            for (int ibatch2 = 0; ibatch2 != numshell2; ++ibatch2) {
              const int offset2 = coffset2[ibatch2]; 
              shared_ptr<Shell> b2 = shell2[ibatch2];

              vector<shared_ptr<Shell> > input;
              input.push_back(b3);
              input.push_back(b2);
              input.push_back(b1);
              input.push_back(b0);

              // pointer to stack
              GradBatch gradbatch(input, 0.0);
              gradbatch.compute();
              const size_t block = gradbatch.size_block();

              // unfortunately the convention is different...
              int jatom[4] = {-1, iatom2, iatom1, iatom0};
              if (gradbatch.swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
              if (gradbatch.swap01()) swap(jatom[0], jatom[1]);
              if (gradbatch.swap23()) swap(jatom[2], jatom[3]);

              for (int i = 0; i != 12; ++i) {
                // if this is a dummy atom
                if (jatom[i/3] < 0) continue;

                const double* ppt = gradbatch.data() + i*block;
                double sum = 0.0;
                for (int j0 = offset0; j0 != offset0 + b0->nbasis(); ++j0) {  
                  for (int j1 = offset1; j1 != offset1 + b1->nbasis(); ++j1) {  
                    for (int j2 = offset2; j2 != offset2 + b2->nbasis(); ++j2, ++ppt) {  
                      sum += *ppt * *o->ptr(j2, j1, j0);
                    }
                  }
                }
                grad[3*jatom[i/3]+i%3] += sum * (iatom0 == iatom1 && ibatch0 == ibatch1 ? 1.0 : 2.0);
              }
            }
          }
        }

      }
    }
  }

  return grad;

}


vector<double> GradEval_base::contract_grad2e_2index(const unique_ptr<double[]>& o) const {
  vector<double> grad(geom_->natom()*3, 0.0);

  // using symmetry (b0 <-> b1)
  vector<shared_ptr<Atom> > aux_atoms = geom_->aux_atoms();
  vector<vector<int> > aux_offsets = geom_->aux_offsets();
  const vector<shared_ptr<Atom> > atoms = geom_->atoms(); 
  const vector<vector<int> > offsets = geom_->offsets();

  for (int iatom0 = 0; iatom0 != geom_->natom(); ++iatom0) {
    const shared_ptr<Atom> catom0 = aux_atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = aux_offsets[iatom0];
    const vector<shared_ptr<Shell> > shell0 = catom0->shells();
    for (int iatom1 = iatom0; iatom1 != geom_->natom(); ++iatom1) {
      const shared_ptr<Atom> catom1 = aux_atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = aux_offsets[iatom1];
      const vector<shared_ptr<Shell> > shell1 = catom1->shells();
      // dummy shell
      const shared_ptr<Shell> b3(new Shell(shell1.front()->spherical()));
      for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
        const int offset0 = coffset0[ibatch0]; 
        shared_ptr<Shell> b0 = shell0[ibatch0];

        for (int ibatch1 = (iatom0 == iatom1 ? ibatch0 : 0); ibatch1 != numshell1; ++ibatch1) {
          const int offset1 = coffset1[ibatch1]; 
          shared_ptr<Shell> b1 = shell1[ibatch1];

          vector<shared_ptr<Shell> > input;
          input.push_back(b1);
          input.push_back(b3);
          input.push_back(b0);
          input.push_back(b3);

          // pointer to stack
          GradBatch gradbatch(input, 0.0);
          gradbatch.compute();
          const size_t block = gradbatch.size_block();

          // unfortunately the convention is different...
          int jatom[4] = {iatom1, -1, iatom0, -1};
          if (gradbatch.swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
          if (gradbatch.swap01()) swap(jatom[0], jatom[1]);
          if (gradbatch.swap23()) swap(jatom[2], jatom[3]);

          for (int i = 0; i != 12; ++i) {
            // if this is a dummy atom
            if (jatom[i/3] < 0) continue;

            const double* ppt = gradbatch.data() + i*block;
            double sum = 0.0;
            for (int j0 = offset0; j0 != offset0 + b0->nbasis(); ++j0) {  
              for (int j1 = offset1; j1 != offset1 + b1->nbasis(); ++j1, ++ppt) {  
                sum += *ppt * o[j1+geom_->naux()*j0];
              }
            }
            grad[3*jatom[i/3]+i%3] -= sum * 0.5 * (iatom0 == iatom1 && ibatch0 == ibatch1 ? 1.0 : 2.0);
          }
        }
      }
    }
  }
  return grad;
}
