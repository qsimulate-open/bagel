//
// BAGEL - Parallel electron correlation program.
// Filename: breitterm.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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
#include <src/rel/breitterm.h>

using namespace std;
using namespace bagel;

BreitTerm::BreitTerm(shared_ptr<const Breit> breit, list<shared_ptr<DFData>> dfdata, list<shared_ptr<const ZMatrix>> cd, vector<int> cd_comp) : breit_(breit) {

  int dat = 0;
  for (auto& i : dfdata) {
    for (auto& j : i->basis()) {
      for (int k = 0; k != breit_->data().size(); ++k) {
        int m = 0;
        for (auto& l : cd) {
          if (breit_->index(k).first == j->comp() && breit_->index(k).second == cd_comp[m]) {
            shared_ptr<ZMatrix> s12a(new ZMatrix(*(i->df()->data2())));
            shared_ptr<ZMatrix> breit_cd(new ZMatrix(*(breit_->data(k))));
            data_[dat].push_back(shared_ptr<ZMatrix> (new ZMatrix(*s12a * *breit_cd * *(l))));
          }
          m++;
        }
      }
    }
    dat++;
  }
#if 0
  shared_ptr<ZMatrix> sum = btdata.front()->clone();
  auto btiter = btdata.begin();
  for (auto& i : dfc) { 
    for (auto& j : i->basis())
      sum->zaxpy(j->fac(), *btiter++);
  }
  assert(btiter == btdata.end());
#endif

}

