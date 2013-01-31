//
// BAGEL - Parallel electron correlation program.
// Filename: dfock.cc
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


#include <src/util/constants.h>
#include <src/rel/dfock.h>
#include <src/util/matrix.h>

using namespace std;
using namespace bagel;

void DFock::two_electron_part(const array<shared_ptr<const ZMatrix>, 4> ocoeff, const bool rhf, const double scale_exchange) {

  if (!rhf) throw logic_error("DFock::two_electron_part() is not implemented for non RHF cases");

  const int n = geom_->nbasis();
  const int nele = geom_->nele();
  shared_ptr<const DFDist> df = geom_->df();
  shared_ptr<const DFDist> dfs_total = geom_->dfs();

  // get individual df dist objects for each block
  vector<shared_ptr<DFDist>> dfs = dfs_total->split_blocks();

  // Separate Coefficients into real and imaginary
  array<shared_ptr<const Matrix>, 4> rocoeff;
  array<shared_ptr<const Matrix>, 4> iocoeff;
  array<shared_ptr<const Matrix>, 4> trocoeff;
  array<shared_ptr<const Matrix>, 4> tiocoeff;

  for (int i = 0; i != 4; ++i) {
    rocoeff[i] = ocoeff[i]->get_real_part();
    iocoeff[i] = ocoeff[i]->get_imag_part();
    trocoeff[i] = rocoeff[i]->transpose();
    tiocoeff[i] = iocoeff[i]->transpose();
  }

  compute_half_complex(rocoeff, iocoeff, df, dfs);

  complex<double> coeffi(0.0, 1.0);
  complex<double> coeff1(1.0, 0.0);

  for (int i = 0; i != 2; ++i) {
    const int j = large_half_[i]->basis().second;
    array<shared_ptr<Matrix>, 2> large_data = large_half_[i]->compute_Jop(df, trocoeff[j], tiocoeff[j], large_half_[i]->basis(), large_half_[i]->coord());
    // J should be added to both blocks
    add_real_block(coeff1, 0, 0, n, n, large_data[0]);
    add_real_block(coeffi, 0, 0, n, n, large_data[1]);
    add_real_block(coeff1, n, n, n, n, large_data[0]);
    add_real_block(coeffi, n, n, n, n, large_data[1]);

    for(int j = i; j != 2; ++j) {
      array<shared_ptr<Matrix>, 2> large_data2 = large_half_[i]->form_2index(large_half_[j]);
      add_real_block(coeff1, n*large_half_[i]->basis().second, n*large_half_[j]->basis().second, n, n, large_data2[0]);
      add_real_block(coeffi, n*large_half_[i]->basis().second, n*large_half_[j]->basis().second, n, n, large_data2[1]);
      if (j != i) {
        add_real_block(coeff1, n*large_half_[j]->basis().second, n*large_half_[i]->basis().second, n, n, large_data2[0]->transpose());
        add_real_block(coeffi, n*large_half_[j]->basis().second, n*large_half_[i]->basis().second, n, n, (*large_data2[1]->transpose() * (-1.0)).data());
      }
    }
  }

  // Small Half Transforms; swapped only needs xy, xz, and yz
  const double tc = 1.0 / (2.0*c__);
  const complex<double> sscoeff1 = coeff1 * ::pow(tc,4.0);
  const complex<double> sscoeffi = coeffi * ::pow(tc,4.0);

#if 1
  for (int i = 0; i != small_half_.size(); ++i) {
    for (int j = i; j != small_half_.size(); ++j) {
      array<shared_ptr<Matrix>, 2> small_data = small_half_[i]->form_2index(small_half_[j]);
      add_real_block(sscoeff1, n*(2+small_half_[i]->basis().second), n*(2+small_half_[j]->basis().second), n, n, small_data[0]);
      add_real_block(sscoeffi, n*(2+small_half_[i]->basis().second), n*(2+small_half_[j]->basis().second), n, n, small_data[1]);
      if (j != i) {
        add_real_block(sscoeff1, n*(2+small_half_[j]->basis().second), n*(2+small_half_[i]->basis().second), n, n, small_data[0]->transpose());
        add_real_block(sscoeffi, n*(2+small_half_[j]->basis().second), n*(2+small_half_[i]->basis().second), n, n, (*small_data[1]->transpose() * (-1.0)).data());
      }
    }
  }
#endif
  
#if 1
  pair<const int, const int> aa = make_pair(0,0);
  pair<const int, const int> ab = make_pair(0,1);
  pair<const int, const int> ba = make_pair(1,0);
  pair<const int, const int> bb = make_pair(1,1);
  pair<const int, const int> l_coord = make_pair(-1, -1);

  pair<const int, const int> xx = make_pair(0,0);
  pair<const int, const int> xy = make_pair(0,1);
  pair<const int, const int> xz = make_pair(0,2);
  pair<const int, const int> yx = make_pair(1,0);
  pair<const int, const int> yy = make_pair(1,1);
  pair<const int, const int> yz = make_pair(1,2);
  pair<const int, const int> zx = make_pair(2,0);
  pair<const int, const int> zy = make_pair(2,1);
  pair<const int, const int> zz = make_pair(2,2);
  pair<const int, const int> eri_coord[6] = {xx, xy, xz, yy, yz, zz};
  pair<const int, const int> eri_basis1[6] = {aa, aa, ab, aa, ab, aa};
  pair<const int, const int> eri_basis2[6] = {bb, bb, ba, bb, ba, bb};

#if 1
  for (int i = 0; i != dfs.size(); ++i) {
    for (int j = 0; j != small_half_.size(); ++j) {
      const int k = small_half_[j]->basis().second+2;
      array<shared_ptr<Matrix>, 2> small_data1 = small_half_[j]->compute_Jop(dfs[i], trocoeff[k], tiocoeff[k], eri_basis1[i], eri_coord[i]); 
      // TODO target1 and target 2 should be obtained by a better way
      const tuple<int,int> target1 = make_tuple(eri_basis1[i].first+2, eri_basis1[i].second+2);
      const tuple<int,int> target2 = make_tuple(eri_basis2[i].first+2, eri_basis2[i].second+2);

      const complex<double> sign = (i == 1 || i == 2) ? -coeff1 : coeff1;

      add_real_block(sscoeff1, n*get<0>(target1), n*get<1>(target1), n, n, small_data1[0]);
      add_real_block(sscoeffi, n*get<0>(target1), n*get<1>(target1), n, n, small_data1[1]);
      add_real_block(sscoeff1*sign, n*get<0>(target2), n*get<1>(target2), n, n, small_data1[0]);
      add_real_block(sscoeffi*sign, n*get<0>(target2), n*get<1>(target2), n, n, small_data1[1]);

      if (i == 1 || i == 2 || i == 4) {
        add_real_block(sscoeff1, n*get<1>(target1), n*get<0>(target1), n, n, (*small_data1[0]->transpose() * (i == 2 ? 1.0 : -1.0)).data());
        add_real_block(sscoeffi, n*get<1>(target1), n*get<0>(target1), n, n, (*small_data1[1]->transpose() * (i == 2 ? 1.0 : -1.0)).data());
        add_real_block(sscoeff1*sign, n*get<1>(target2), n*get<0>(target2), n, n, (*small_data1[0]->transpose() * (i == 2 ? 1.0 : -1.0)).data());
        add_real_block(sscoeffi*sign, n*get<1>(target2), n*get<0>(target2), n, n, (*small_data1[1]->transpose() * (i == 2 ? 1.0 : -1.0)).data());
      }
    }
  }
#endif
#endif

#if 0
  // Small - Large Transforms
  for (int i = 0; i != 2; ++i) {
    for (int j = 0; j != small_half_.size(); ++j) {
      array<shared_ptr<Matrix>, 2> large_small_data;
      large_small_data = large_half[i]->form_2index_small(small_half_[j]);
      //TODO figure out what to do with this damn data
    }
  }
#endif

}

void DFock::compute_half_complex(array<shared_ptr<const Matrix>, 4> rocoeff, array<shared_ptr<const Matrix>, 4> iocoeff, shared_ptr<const DFDist> df, 
                                 vector<shared_ptr<DFDist> > dfs) {
  const int a_basis = 2;
  const int b_basis = 3;
  pair<const int, const int> aa = make_pair(0,0);
  pair<const int, const int> ab = make_pair(0,1);
  pair<const int, const int> ba = make_pair(1,0);
  pair<const int, const int> bb = make_pair(1,1);
  pair<const int, const int> l_coord = make_pair(-1, -1);

  pair<const int, const int> xx = make_pair(0,0);
  pair<const int, const int> xy = make_pair(0,1);
  pair<const int, const int> xz = make_pair(0,2);
  pair<const int, const int> yx = make_pair(1,0);
  pair<const int, const int> yy = make_pair(1,1);
  pair<const int, const int> yz = make_pair(1,2);
  pair<const int, const int> zx = make_pair(2,0);
  pair<const int, const int> zy = make_pair(2,1);
  pair<const int, const int> zz = make_pair(2,2);

  large_half_[0] = shared_ptr<DFHalfComplex>(new DFHalfComplex(df, rocoeff[0], iocoeff[0], false, l_coord, aa));
  large_half_[1] = shared_ptr<DFHalfComplex>(new DFHalfComplex(df, rocoeff[1], iocoeff[1], false, l_coord, bb));

  // XX
  small_half_[0] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[0], rocoeff[a_basis], iocoeff[a_basis], false, xx, aa));
  small_half_[1] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[0], rocoeff[b_basis], iocoeff[b_basis], false, xx, bb));

  // XY
  small_half_[2] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[a_basis], iocoeff[a_basis], false, xy, aa));
  small_half_[3] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[b_basis], iocoeff[b_basis], false, xy, bb));

  // XZ
  small_half_[4] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[a_basis], iocoeff[a_basis], false, xz, ab));
  small_half_[5] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[b_basis], iocoeff[b_basis], false, xz, ba));

  // YX
  small_half_[6] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[a_basis], iocoeff[a_basis], true, yx, aa));
  small_half_[7] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[b_basis], iocoeff[b_basis], true, yx, bb));

  // YY
  small_half_[8] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[3], rocoeff[a_basis], iocoeff[a_basis], false, yy, aa));
  small_half_[9] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[3], rocoeff[b_basis], iocoeff[b_basis], false, yy, bb));

  // YZ
  small_half_[10] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[a_basis], iocoeff[a_basis], false, yz, ab));
  small_half_[11] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[b_basis], iocoeff[b_basis], false, yz, ba));

  // ZX
  small_half_[12] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[a_basis], iocoeff[a_basis], true, zx, ab));
  small_half_[13] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[b_basis], iocoeff[b_basis], true, zx, ba));

  // ZY
  small_half_[14] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[a_basis], iocoeff[a_basis], true, zy, ab));
  small_half_[15] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[b_basis], iocoeff[b_basis], true, zy, ba));

  // ZZ
  small_half_[16] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[5], rocoeff[a_basis], iocoeff[a_basis], false, zz, aa));
  small_half_[17] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[5], rocoeff[b_basis], iocoeff[b_basis], false, zz, bb));
}

