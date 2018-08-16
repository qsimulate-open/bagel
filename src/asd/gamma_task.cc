//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gamma_task.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/asd/gamma_task.h>

using namespace std;
using namespace bagel;


//TODO template to other VecType
template <>
void GammaTask<CASDvec>::compute() {
  constexpr int nops = 4;
  const int norb = tree_->norb();

  auto action = [] (const int op) { return is_creation(GammaSQ(op)); };
  auto spin = [] (const int op) { return is_alpha(GammaSQ(op)); };

  shared_ptr<GammaBranch<CASDvec>> first = tree_->base()->branch(operation_);
  assert(first->active()); // This should have been checked before sending it to the TaskQueue

  shared_ptr<CASDvec> avec = tree_->ket()->apply_and_allocate(action(static_cast<int>(operation_)), spin(static_cast<int>(operation_))); //workspace
  avec->apply_and_fill(tree_->ket(), a_, action(static_cast<int>(operation_)), spin(static_cast<int>(operation_))); //fill

  for (auto& ibra : first->bras())
    dot_product(ibra.second, avec, first->gammas().find(ibra.first)->second->element_ptr(0,a_));

  for (int j = 0; j < nops; ++j) {
    auto second = first->branch(j);
    if (!second->active()) continue;

    shared_ptr<CASDvec> bvec = avec->apply_and_allocate(action(j), spin(j));

    for (int b = 0; b < norb; ++b) {
      if (b==a_ && j==static_cast<int>(operation_)) continue;
      bvec->apply_and_fill(avec, b, action(j), spin(j));
      for (auto& jbra : second->bras())
        dot_product(jbra.second, bvec, second->gammas().find(jbra.first)->second->element_ptr(0, a_*norb + b));

      for (int k = 0; k < nops; ++k) {
        shared_ptr<GammaBranch<CASDvec>> third = second->branch(k);
        if (!third->active()) continue;

        shared_ptr<CASDvec> cvec = bvec->apply_and_allocate(action(k), spin(k));

        for (int c = 0; c < norb; ++c) {
          if (b==c && k==j) continue;
          cvec->apply_and_fill(bvec, c, action(k), spin(k));
          for (auto& kbra : third->bras())
            dot_product(kbra.second, cvec, third->gammas().find(kbra.first)->second->element_ptr(0, a_*norb*norb + b*norb + c));
        }
      }
    }
  }

}
