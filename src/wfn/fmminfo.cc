//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fmminfo.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#include <src/wfn/fmminfo.h>

using namespace std;
using namespace bagel;

FMMInfo::FMMInfo(const vector<shared_ptr<const Atom>>& atoms, const vector<vector<int>>& offsets, const string type) {

  const string extent_type = (type == "") ? "yang" : type;

  vector<int> wk;
  vector<shared_ptr<const Shell>> basis;
  for (int n = 0; n != atoms.size(); ++n) {
    const vector<int> tmpoff = offsets[n];
    wk.insert(wk.end(), tmpoff.begin(), tmpoff.end());
    const vector<shared_ptr<const Shell>> tmpsh = atoms[n]->shells();
    basis.insert(basis.end(), tmpsh.begin(), tmpsh.end());
  }
  const int nsh = basis.size();

  shellpairs_.resize(nsh * nsh);
  for (int i0 = 0; i0 != nsh; ++i0) {
    for (int i1 = 0; i1 != nsh; ++i1) {
      const int i01 = i0 * nsh + i1;
      shellpairs_[i01] = make_shared<const ShellPair>(array<shared_ptr<const Shell>, 2>{{basis[i1], basis[i0]}},
                                                      array<int, 2>{{wk[i1], wk[i0]}}, make_pair(i1, i0), extent_type);
    }
  }
}
