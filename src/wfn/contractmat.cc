//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: contractmat.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#include <src/wfn/contractmat.h>
#include <src/integral/os/overlapbatch.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(ContractMat)


ContractMat::ContractMat(shared_ptr<const Molecule> mol, int nunc) : Matrix(nunc, mol->nbasis()) {
  init(mol);
}


void ContractMat::init(shared_ptr<const Molecule> mol) {
  int npos = 0, mpos = 0, nsize = 0, msize = 0;
  auto oa = mol->offsets().begin();
  for (auto a = mol->atoms().begin(); a != mol->atoms().end(); ++a, ++oa) {
    auto ob = oa->begin();
    for (auto b = (*a)->shells().begin(); b != (*a)->shells().end(); ++b, ++ob) {
      vector<vector<double>> contractions((*b)->contractions().size());
      auto citer = contractions.begin();
      for (auto iter = (*b)->contractions().begin(); iter != (*b)->contractions().end(); ++iter, ++citer) {
        auto eiter = (*b)->exponents().begin();
        double denom = 1.0;
        for (int ii = 2; ii <= (*b)->angular_number(); ++ii) denom *= 2 * ii - 1;
        for (auto diter = iter->begin(); diter != iter->end(); ++diter, ++eiter) {
          citer->push_back(*diter * std::sqrt(denom) / (pow(2.0 * *eiter / pi__, 0.75) * pow(std::sqrt(4.0 * *eiter), (*b)->angular_number())));
        }
      }
      int duplicates = (*b)->spherical() ? 2 * (*b)->angular_number() + 1 : ((*b)->angular_number() + 1) * ((*b)->angular_number() + 2) / 2;
      const vector<pair<int, int>> contraction_ranges = (*b)->contraction_ranges();
      for (int r = 0; r != contraction_ranges.size(); r += msize) {
        msize = 0;
        for (int i = r; i != contraction_ranges.size() && contraction_ranges[i].first == contraction_ranges[r].first; ++i) ++msize;
        nsize = contraction_ranges[r].second - contraction_ranges[r].first;
        for (int i = 0; i != msize; ++i) {
          for (int j = 0; j != nsize; ++j) {
            for (int a = 0; a != duplicates; ++a) {
              assert(npos + j * duplicates + a < ndim() && mpos + i * duplicates + a < mdim());
              element(npos + j * duplicates + a, mpos + i * duplicates + a) = contractions[r + i][contraction_ranges[r].first + j];
            }
          }
        }
        npos += nsize * duplicates;
        mpos += msize * duplicates;
      }
    }
  }
}
