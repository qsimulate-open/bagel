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
#include <src/util/taskqueue.h>
#include <src/util/resources.h>
#include <array>

using namespace std;


shared_ptr<GradFile> GradEval_base::contract_gradient(const shared_ptr<const Matrix1e> d, const shared_ptr<const Matrix1e> w,
                                                      const shared_ptr<const DF_AO> o, const unique_ptr<double[]>& o2) {

  shared_ptr<Grad1eFile> g1(new Grad1eFile(geom_));
  shared_ptr<Grad1eFile> go(new Grad1eFile(geom_));
  compute_grad1e_integrals(g1, go);

  vector<double> grad = contract_grad1e(d, w, g1, go); 
  vector<GradTask> tasks = contract_grad2e(o);
  vector<GradTask> tmp = contract_grad2e_2index(o2);
  tasks.insert(tasks.end(), tmp.begin(), tmp.end());

  TaskQueue<GradTask> tq(tasks);
  tq.compute(resources__->max_num_threads());

  shared_ptr<GradFile> g0(new GradFile(grad));
  return shared_ptr<GradFile>(new GradFile(*g0+*grad_));
}


void GradEval_base::compute_grad1e_integrals(shared_ptr<Grad1eFile> g1, shared_ptr<Grad1eFile> go) const {

  const int natom = geom_->natom();

  const vector<shared_ptr<const Atom> > atoms = geom_->atoms(); 
  const vector<vector<int> > offsets = geom_->offsets();
  const int nbasis = geom_->nbasis();

  // TODO perhaps we could reduce operation by a factor of 2
  for (int iatom0 = 0; iatom0 != natom; ++iatom0) {
    const shared_ptr<const Atom> catom0 = atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = offsets[iatom0];
    const vector<shared_ptr<const Shell> > shell0 = catom0->shells();

    for (int iatom1 = 0; iatom1 != natom; ++iatom1) {
      const shared_ptr<const Atom> catom1 = atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = offsets[iatom1];
      const vector<shared_ptr<const Shell> > shell1 = catom1->shells();

      for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
        const int offset0 = coffset0[ibatch0]; 
        shared_ptr<const Shell> b0 = shell0[ibatch0];
        for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
          const int offset1 = coffset1[ibatch1]; 
          shared_ptr<const Shell> b1 = shell1[ibatch1];
          vector<shared_ptr<const Shell> > input = {{b1, b0}};

          const int dimb1 = input[0]->nbasis(); 
          const int dimb0 = input[1]->nbasis(); 

          {
            GNAIBatch batch2(input, geom_, tie(iatom1, iatom0));
            batch2.compute();
            const double* ndata = batch2.data();
            const size_t s = batch2.size_block();
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
            const size_t s = batch.size_block();
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
            const size_t s = batch.size_block();
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




vector<GradTask> GradEval_base::contract_grad2e(const shared_ptr<const DF_AO> o) {
  vector<GradTask> out;

  // loop over atoms (using symmetry b0 <-> b1)
  int iatom0 = 0;
  auto oa0 = geom_->offsets().begin();
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++oa0, ++iatom0) {
    int iatom1 = iatom0;
    auto oa1 = oa0;
    for (auto a1 = a0; a1 != geom_->atoms().end(); ++a1, ++oa1, ++iatom1) {
      int iatom2 = 0;
      auto oa2 = geom_->aux_offsets().begin(); 
      for (auto a2 = geom_->aux_atoms().begin(); a2 != geom_->aux_atoms().end(); ++a2, ++oa2, ++iatom2) {

        // dummy shell
        const shared_ptr<const Shell> b3(new Shell((*a0)->shells().front()->spherical()));

        auto o0 = oa0->begin();
        for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
          auto o1 = a0!=a1 ? oa1->begin() : o0;
          for (auto b1 = (a0!=a1 ? (*a1)->shells().begin() : b0); b1 != (*a1)->shells().end(); ++b1, ++o1) {
            auto o2 = oa2->begin();
            for (auto b2 = (*a2)->shells().begin(); b2 != (*a2)->shells().end(); ++b2, ++o2) {
              vector<shared_ptr<const Shell> > input = {{b3, *b2, *b1, *b0}};
              vector<int> atoms = {{iatom0, iatom1, iatom2}};
              vector<int> offs = {{*o0, *o1, *o2}};

              GradTask task(input, atoms, offs, o, this); 
              out.push_back(task);
            }
          }
        }

      }
    }
  }
  return out;
}


vector<GradTask> GradEval_base::contract_grad2e_2index(const unique_ptr<double[]>& o) {
  vector<GradTask> out;

  shared_ptr<Geometry> auxgeom(new Geometry(geom_->aux_atoms(), multimap<string,string>()));
  shared_ptr<Matrix1e> den(new Matrix1e(auxgeom));
  copy(o.get(), o.get()+geom_->naux()*geom_->naux(), den->data());

  // using symmetry (b0 <-> b1)
  int iatom0 = 0;
  auto oa0 = geom_->aux_offsets().begin();
  for (auto a0 = geom_->aux_atoms().begin(); a0 != geom_->aux_atoms().end(); ++a0, ++oa0, ++iatom0) {
    int iatom1 = iatom0;
    auto oa1 = oa0;
    for (auto a1 = a0; a1 != geom_->aux_atoms().end(); ++a1, ++oa1, ++iatom1) {

      // dummy shell
      const shared_ptr<const Shell> b3(new Shell((*a0)->shells().front()->spherical()));

      auto o0 = oa0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
        auto o1 = a0!=a1 ? oa1->begin() : o0;
        for (auto b1 = (a0!=a1 ? (*a1)->shells().begin() : b0); b1 != (*a1)->shells().end(); ++b1, ++o1) {


          vector<shared_ptr<const Shell> > input = {{*b1, b3, *b0, b3}};
          vector<int> atoms = {{iatom0, iatom1}};
          vector<int> offs = {{*o0, *o1}};

          GradTask task(input, atoms, offs, den, this);
          out.push_back(task);
        }
      }
    }
  }
  return out;
}
