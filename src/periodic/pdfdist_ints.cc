//
// BAGEL - Parallel electron correlation program.
// Filename: pdfdist_ints.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

#include <src/periodic/pdfdist_ints.h>
#include <src/periodic/pdfinttask.h>

using namespace bagel;
using namespace std;

void PDFDist_ints::compute_3index(const vector<shared_ptr<const Shell>>& ashell,
                                  const vector<shared_ptr<const Shell>>& b0shell,
                                  const vector<shared_ptr<const Shell>>& bgshell) {
  Timer time;

  TaskQueue<PDFIntTask> tasks(b0shell.size() * bgshell.size() * ashell.size());
  auto i3 = make_shared<const Shell>(ashell.front()->spherical());

  int j2 = 0;
  for (auto& i2 : bgshell) {
    int j1 = 0;
    for (auto& i1 : b0shell) {
      int j0 = 0;
      for (auto& i0 : ashell) {
        tasks.emplace_back((array<shared_ptr<const Shell>, 4>{{i3, i0, i1, i2}}), (array<int, 3>{{j2, j1, j0}}), block_[0]);
        j0 += i0->nbasis();
      }
      j1 += i1->nbasis();
    }
    j2 += i2->nbasis();
  }

  time.tick_print("3-index integrals prep");
  tasks.compute();
  time.tick_print("3-index integrals");
}
