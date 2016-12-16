//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: jacobi_pairs.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __BAGEL_MATH_JACOBI_PAIRS_H
#define __BAGEL_MATH_JACOBI_PAIRS_H

#include <string>
#include <algorithm>
#include <vector>

namespace bagel {

/************************************************************************************
* Methods to generate Jacobi pairs. See reference:                                  *
*  B. B. Zhou and R. P. Brent, Proc. Euromicro Workshop on Parallel and             *
*     Distributed Processing (San Remo, Italy), IEEE CS Press, 1995, 401-408.       *
************************************************************************************/

class JacobiPairs {
  protected:
    std::vector<std::vector<std::pair<int, int>>> fullsweep_;

  public:
    std::vector<std::vector<std::pair<int, int>>>::iterator begin() { return fullsweep_.begin(); }
    std::vector<std::vector<std::pair<int, int>>>::iterator end() { return fullsweep_.end(); }
};

class JacobiRoundRobin : public JacobiPairs {
  public:
    JacobiRoundRobin(const int nstart, const int nfence) {
      const int norb = nfence - nstart;

      const int ntotal = norb + (norb%2);
      const int nhalf = ntotal/2;

      std::vector<int> orbitals(ntotal, 0);
      for (int i = 0; i < nhalf; ++i) { orbitals[i] = 2*i + nstart; orbitals[i + nhalf] = 2*i + nstart + 1; }
      orbitals.back() *= 1 - 2*(norb%2);

      const int nsubsweeps = ntotal - 1;

      int count = 0;
      for (int isub = 0; isub < nsubsweeps; ++isub) {
        std::vector<std::pair<int,int>> pairlist;
        for (int jpair = 0; jpair < nhalf; ++jpair) {
          const int ii = orbitals[jpair];
          const int jj = orbitals[ntotal - jpair - 1];
          if ( ii < 0 || jj < 0 ) continue;
          pairlist.emplace_back(ii, jj);
          ++count;
        }
        fullsweep_.push_back(pairlist);
        std::rotate(orbitals.begin(), orbitals.begin() + 1, orbitals.end() - 1);
      }
      assert( count == ((norb*(norb-1))/2) );
    }
};

class JacobiOddEven : public JacobiPairs {
  public:
    JacobiOddEven(const int nstart, const int nfence) {
      const int norb = nfence - nstart;
      std::vector<int> orbitals(norb, 0);
      std::iota(orbitals.begin(), orbitals.end(), 0);

      const int nsubsweeps = norb;

      int count = 0;
      for (int isub = 0; isub < nsubsweeps; ++isub) {
        std::vector<std::pair<int, int>> pairlist;
        for (int i = isub%2; i < norb; i+=2) {
          if ( (i+1) < norb ) {
            pairlist.emplace_back(orbitals[i], orbitals[i+1]);
            std::swap(orbitals[i], orbitals[i+1]);
            ++count;
          }
        }
        fullsweep_.push_back(pairlist);
      }

      assert( count == ((norb*(norb-1))/2) );
    }
};

class JacobiRing : public JacobiPairs {
  public:
    JacobiRing(const int nstart, const int nfence) {
      const int norb = nfence - nstart;
      const int ntotal = norb + (norb%2);
      const int nhalf = ntotal/2;

      std::vector<int> upper(nhalf, 0);
      std::vector<int> lower(nhalf, 0);
      for (int i = 0; i < nhalf; ++i) { lower[i] = 2*i + nstart; upper[i] = 2*i + 1 + nstart; }
      upper.back() *= 1 - 2*(norb%2);

      const int nsubsweeps = ntotal - 1;

      int count = 0;
      for (int isub = 0; isub < nsubsweeps; ++isub) {
        std::vector<std::pair<int, int>> pairlist;
        for (int i = 0; i < nhalf; ++i) {
          const int ii = upper[i];
          const int jj = lower[i];
          if ( ii < 0 || jj < 0 ) continue;
          pairlist.emplace_back(ii, jj);
          ++count;
        }
        fullsweep_.push_back(pairlist);

        const int switchat = nhalf - (isub/2 + 1);
        std::swap(upper[switchat], lower[switchat]);
        std::rotate(upper.begin(), upper.begin() + 1, upper.end());
      }
      assert( count == (norb*(norb-1))/2 );
    }
};

}

#endif
