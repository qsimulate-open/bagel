//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: slater2e.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/scf/hf/fock.h>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <src/integral/rys/eribatch.h>
#include <src/integral/rys/slaterbatch.h>
#include <src/util/constants.h>

using namespace std;
using namespace bagel;

typedef std::shared_ptr<Matrix1e> RefDensity;
typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<Shell> RefShell;

#if 0
void Fock::slater_two_electron_part() {
  // for debug
  density_->fill_upper();

  const vector<RefAtom> atoms = geom_->atoms();
  vector<RefShell> basis;
  vector<int> offset;
  int cnt = 0;
  for (vector<RefAtom>::const_iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
    const vector<RefShell> tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());
    const vector<int> tmpoff = geom_->offset(cnt);
    offset.insert(offset.end(), tmpoff.begin(), tmpoff.end());
  }

  const int shift = sizeof(int) * 4;
  const int size = basis.size();

  // first make max_density_change vector for each batch pair.
  const double* density_data = density_->data();

  vector<double> max_density_change(size * size);
  for (int i = 0; i != size; ++i) {
    const int ioffset = offset[i];
    const int isize = basis[i]->nbasis();
    for (int j = i; j != size; ++j) {
      const int joffset = offset[j];
      const int jsize = basis[j]->nbasis();

      double cmax = 0.0;
      for (int ii = ioffset; ii != ioffset + isize; ++ii) {
        const int iin = ii * nbasis_;
        for (int jj = joffset; jj != joffset + jsize; ++jj) {
          cmax = max(cmax, fabs(density_data[iin + jj]));
        }
      }
      const int ij = i * size + j;
      const int ji = i * size + j;
      max_density_change[ij] = cmax;
      max_density_change[ji] = cmax;
    }
  }

  RefPetite plist = geom_->plist();;
  const bool c1 = plist->nirrep() == 1;

  for (int i0 = 0; i0 != size; ++i0) {
    if (!plist->in_p1(i0)) continue;

    const RefShell b0 = basis[i0];
    const int b0offset = offset[i0];
    const int b0size = b0->nbasis();
    for (int i1 = i0; i1 != size; ++i1) {
      const unsigned int i01 = i0 *size + i1;
      if (!plist->in_p2(i01)) continue;

      const RefShell b1 = basis[i1];
      const int b1offset = offset[i1];
      const int b1size = b1->nbasis();

      const double density_change_01 = max_density_change[i01] * 4.0;

      for (int i2 = i0; i2 != size; ++i2) {
        const RefShell b2 = basis[i2];
        const int b2offset = offset[i2];
        const int b2size = b2->nbasis();

        const double density_change_02 = max_density_change[i0 * size + i2];
        const double density_change_12 = max_density_change[i1 * size + i2];

        for (int i3 = i2; i3 != size; ++i3) {
          const unsigned int i23 = i2 * size + i3;
          if (i23 < i01) continue;
          int ijkl = plist->in_p4(i01, i23, i0, i1, i2, i3);
          if (ijkl == 0) continue;

          const double density_change_23 = max_density_change[i2 * size + i3] * 4.0;
          const double density_change_03 = max_density_change[i0 * size + i2];
          const double density_change_13 = max_density_change[i0 * size + i2];

          const bool eqli01i23 = (i01 == i23);

          const RefShell b3 = basis[i3];
          const int b3offset = offset[i3];
          const int b3size = b3->nbasis();

// schwarz prescreening >>>>
          const double mulfactor = max(max(max(density_change_01, density_change_02),
                                           max(density_change_12, density_change_23)),
                                           max(density_change_03, density_change_13));
          const double integral_bound = mulfactor * schwarz_[i01] * schwarz_[i23];
          const bool skip_schwarz = integral_bound < schwarz_thresh__;
#if 0
          if (skip_schwarz) continue;
#endif

          vector<RefShell> input;
          input.push_back(b3);
          input.push_back(b2);
          input.push_back(b1);
          input.push_back(b0);

//#define COMPUTE_ERI
#ifdef COMPUTE_ERI
          ERIBatch eribatch(input, 1.0e100);
          eribatch.compute();
          const double* eridata = eribatch.data();
          assert(!eribatch.data2_exists()) {
          assert((int)eribatch.data_size() == b0size * b1size * b2size * b3size);
#if 0
if (i0 == 1 && i1 == 2 && i2 == 5 && i3 == 6) {
  int i = 0;
  for (int z0 = 0; z0 != b0size; ++z0) {
    for (int z1 = 0; z1 != b1size; ++z1) {
      for (int z2 = 0; z2 != b2size; ++z2) {
        for (int z3 = 0; z3 != b3size; ++z3, ++i) {
          cout << z0 << " " << z1 << " " << z2 << " " << z3 << " " << setprecision(10) << eridata[i] << endl;
        }
      }
    }
  }
  assert(false);
}
#endif
#else
          SlaterBatch slaterbatch(input, 0.0, 1.5, true);
          slaterbatch.compute();
          const double* eridata = slaterbatch.data();
          assert(slaterbatch.data2_exists());
          const double* eridatay = slaterbatch.data2();
#if 0
if (i0 == 1 && i1 == 2 && i2 == 5 && i3 == 6) {
  int i = 0;
  for (int z0 = 0; z0 != b0size; ++z0) {
    for (int z1 = 0; z1 != b1size; ++z1) {
      for (int z2 = 0; z2 != b2size; ++z2) {
        for (int z3 = 0; z3 != b3size; ++z3, ++i) {
          cout << z0 << " " << z1 << " " << z2 << " " << z3 << " " << setprecision(10) << eridatay[i] << endl;
        }
      }
    }
  }
  assert(false);
}
#endif
#endif

#if 1

          for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
            const int j0n = j0 * nbasis_;

            for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
              const unsigned int nj01 = (j0 << shift) + j1;
              const bool skipj0j1 = (j0 > j1);
              if (skipj0j1) {
                eridata += b2size * b3size;
                continue;
              }

              const bool eqlj0j1 = (j0 == j1);
              const double scal01 = (eqlj0j1 ? 0.5 : 1.0) * static_cast<double>(ijkl);
              const int j1n = j1 * nbasis_;

              for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                const int maxj1j2 = max(j1, j2);
                const int minj1j2 = min(j1, j2);
                const int minj1j2n = minj1j2 * nbasis_;

                const int maxj0j2 = max(j0, j2);
                const int minj0j2 = min(j0, j2);
                const int minj0j2n = minj0j2 * nbasis_;
                const int j2n = j2 * nbasis_;

                for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                  const bool skipj2j3 = (j2 > j3);
                  const unsigned int nj23 = (j2 << shift) + j3;
                  const bool skipj01j23 = (nj01 > nj23) && eqli01i23;

                  if (skipj2j3 || skipj01j23) continue;

                  const int maxj1j3 = max(j1, j3);
                  const int minj1j3 = min(j1, j3);

                  double intval = *eridata * scal01 * (j2 == j3 ? 0.5 : 1.0) * (nj01 == nj23 ? 0.5 : 1.0);
                  const double intval4 = 4.0 * intval;

                  data_[j0n + j1] += density_data[j2n + j3] * intval4;
                  data_[j2n + j3] += density_data[j0n + j1] * intval4;
                  data_[j0n + j3] -= density_data[j1n + j2] * intval;
                  data_[minj1j2n + maxj1j2] -= density_data[j0n + j3] * intval;
                  data_[minj0j2n + maxj0j2] -= density_data[j1n + j3] * intval;
                  data_[minj1j3 * nbasis_ + maxj1j3] -= density_data[j0n + j2] * intval;
                }
              }
            }
          }
#endif

        }
      }
    }
  }
}

#endif
