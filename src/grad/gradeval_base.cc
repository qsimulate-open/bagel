//
// BAGEL - Parallel electron correlation program.
// Filename: gradeval_base.cc
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


#include <src/grad/gradeval_base.h>
#include <src/util/taskqueue.h>
#include <src/util/resources.h>
#include <array>

using namespace std;
using namespace bagel;

shared_ptr<GradFile> GradEval_base::contract_gradient(const shared_ptr<const Matrix> d, const shared_ptr<const Matrix> w,
                                                      const shared_ptr<const DFDist> o, const shared_ptr<const Matrix> o2) {

  vector<GradTask> task  = contract_grad1e(d, w);
  vector<GradTask> task2 = contract_grad2e(o);
  vector<GradTask> task3 = contract_grad2e_2index(o2);
  task.insert(task.end(), task2.begin(), task2.end());
  task.insert(task.end(), task3.begin(), task3.end());

  TaskQueue<GradTask> tq(task);
  tq.compute(resources__->max_num_threads());

  vector<double> gnuc = geom_->compute_grad_vnuc();
  GradFile nuc(gnuc);
  return shared_ptr<GradFile>(new GradFile(*grad_+nuc));
}


vector<GradTask> GradEval_base::contract_grad1e(const shared_ptr<const Matrix> d, const shared_ptr<const Matrix> w) {
  vector<GradTask> out;
  size_t nshell = 0;
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0) nshell += (*a0)->shells().size();
  out.reserve(nshell*nshell);

  // TODO perhaps we could reduce operation by a factor of 2
  int iatom0 = 0;
  auto oa0 = geom_->offsets().begin();
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++oa0, ++iatom0) {
    int iatom1 = 0;
    auto oa1 = geom_->offsets().begin();
    for (auto a1 = geom_->atoms().begin(); a1 != geom_->atoms().end(); ++a1, ++oa1, ++iatom1) {

      auto o0 = oa0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
        auto o1 = oa1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++o1) {

          array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          vector<int> atom = {iatom0, iatom1};
          vector<int> offset_ = {*o0, *o1};

          GradTask task(input, atom, offset_, d, w, this);
          out.push_back(task);
        }
      }
    }
  }
  return out;
}


vector<GradTask> GradEval_base::contract_grad2e(const shared_ptr<const DFDist> o) {
  vector<GradTask> out;
  size_t nshell = 0;
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0) nshell += (*a0)->shells().size();
  size_t nshell2 = 0;
  for (auto a0 = geom_->aux_atoms().begin(); a0 != geom_->aux_atoms().end(); ++a0) nshell2 += (*a0)->shells().size();
  out.reserve(nshell*(nshell+1)*nshell2/2);

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
              array<shared_ptr<const Shell>,4> input = {{b3, *b2, *b1, *b0}};
              vector<int> atoms = {iatom0, iatom1, iatom2};
              vector<int> offs = {*o0, *o1, *o2};

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


vector<GradTask> GradEval_base::contract_grad2e_2index(const shared_ptr<const Matrix> den) {
  vector<GradTask> out;
  size_t nshell2 = 0;
  for (auto a0 = geom_->aux_atoms().begin(); a0 != geom_->aux_atoms().end(); ++a0) nshell2 += (*a0)->shells().size();
  out.reserve(nshell2*(nshell2+1)/2);

  shared_ptr<Geometry> auxgeom(new Geometry(geom_->aux_atoms(), multimap<string,string>()));

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


          array<shared_ptr<const Shell>,4> input = {{*b1, b3, *b0, b3}};
          vector<int> atoms = {iatom0, iatom1};
          vector<int> offs = {*o0, *o1};

          GradTask task(input, atoms, offs, den, this);
          out.push_back(task);
        }
      }
    }
  }
  return out;
}
