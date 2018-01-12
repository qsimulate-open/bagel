//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: contrcoeff.cc
// Copyright (C) 2016 Toru Shiozaki
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


#include <src/mat1e/contrcoeff.h>
#include <src/integral/os/overlapbatch.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(ContrCoeff)


ContrCoeff::ContrCoeff(shared_ptr<const Molecule> mol, int nunc) : Matrix(nunc, mol->nbasis()) {
  init(mol);
}


void ContrCoeff::init(shared_ptr<const Molecule> mol) {
  // int npos = 0, mpos = 0 nsize = 0, msize = 0;
  // auto oa = mol->offsets().begin();
  // for (auto a = mol->atoms().begin(); a != mol->atoms().end(); ++a, ++oa) {
  //   auto ob = oa->begin();
  //   for (auto b = (*a)->shells().begin(); b != (*a)->shells().end(); ++b, ++ob) {
  //     const vector<vector<double>> contractions = (*b)->contractions();
  //     msize = contractions.size();
  //     nsize = contractions[msize - 1].size();
  //     for (int a = 0; a != 2 * (*b)->angular_number() + 1; ++a) {
  //       for (int i = 0; i != msize; ++i) {
  //         for (int j = 0; j != nsize; ++j) {
  //           assert(npos + j < ndim() && mpos + i < mdim());
  //           element(npos + j, mpos + i) = contractions[i][j];
  //         }
  //       }
  //       npos += nsize;
  //       mpos += msize;
  //     }
  //   }
  // }

  int npos = 0, mpos = 0, nsize = 0, msize = 0;
  auto oa = mol->offsets().begin();
  for (auto a = mol->atoms().begin(); a != mol->atoms().end(); ++a, ++oa) {
    auto ob = oa->begin();
    for (auto b = (*a)->shells().begin(); b != (*a)->shells().end(); ++b, ++ob) {

      // const vector<vector<double>> contractions = (*b)->contractions();
      vector<vector<double>> contractions((*b)->contractions().size());
      auto uiter = contractions.begin();
      const vector<pair<int, int>> contraction_ranges = (*b)->contraction_ranges();
      auto citer = (*b)->contraction_ranges().begin();
      for (auto iter = (*b)->contractions().begin(); iter != (*b)->contractions().end(); ++iter, ++uiter, ++citer) {

        auto eiter = (*b)->exponents().begin();
        double denom = 1.0;
        for (int ii = 2; ii <= (*b)->angular_number(); ++ii) denom *= 2 * ii - 1;
        for (auto diter = iter->begin(); diter != iter->end(); ++diter, ++eiter) {
          uiter->push_back(*diter * std::sqrt(denom) / (pow(2.0 * *eiter / pi__, 0.75) * pow(std::sqrt(4.0 * *eiter), (*b)->angular_number())));
        }

        // if ((*a)->basis() != "molden") {
        //   vector<vector<double>> cont {*iter};
        //   vector<pair<int, int>> cran {*citer};
        //   auto current = make_shared<const Shell>((*b)->spherical(), (*b)->position(), (*b)->angular_number(), (*b)->exponents(), cont, cran);
        //   array<shared_ptr<const Shell>,2> cinp {{ current, current }};
        //   OverlapBatch coverlap(cinp);
        //   coverlap.compute();
        //   const double scal = 1.0 / std::sqrt((coverlap.data())[0]);
        //   for (auto& d : *iter) {
        //     cout << "d " << d << " scal " << scal << endl;
        //     uiter->push_back(d / scal);
        //   }
        // }
      }

      for (int r = 0; r != contraction_ranges.size(); r += msize) {
        msize = 0;
        for (int i = r; i != contraction_ranges.size() && contraction_ranges[i].first == contraction_ranges[r].first; ++i) ++msize;
        nsize = contraction_ranges[r].second - contraction_ranges[r].first;
        for (int a = 0; a != 2 * (*b)->angular_number() + 1; ++a) {
          for (int i = 0; i != msize; ++i) {
            for (int j = 0; j != nsize; ++j) {
              assert(npos + j < ndim() && mpos + i < mdim());
              element(npos + j, mpos + i) = contractions[r + i][contraction_ranges[r].first + j];
            }
          }
          npos += nsize;
          mpos += msize;
        }
      }
    }
  }
}
