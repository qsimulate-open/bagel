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
#include <src/rel/smalleribatch.h>

using namespace std;
using namespace bagel;

void DFock::two_electron_part(const array<shared_ptr<const ZMatrix>, 4> ocoeff, const bool rhf, const double scale_exchange) {

  if (!rhf) throw logic_error("DFock::two_electron_part() is not implemented for non RHF cases");

  complex<double> imag (0.0,1.0);
  shared_ptr<const DFDist> df = geom_->df();
  shared_ptr<const DFDist> dfs_total = geom_->form_fit<DFDist_ints<SmallERIBatch>>(1.0e-8, false); // TODO thresh should be controlled from the input deck

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

  // Large Half Transforms
  array<shared_ptr<DFHalfComplex>, 2> large_half;
  pair<const int, const int> aa = make_pair(0,0);
  pair<const int, const int> ab = make_pair(0,1);
  pair<const int, const int> ba = make_pair(1,0);
  pair<const int, const int> bb = make_pair(1,1);
  pair<const int, const int> l_coord = make_pair(-1, -1);

  large_half[0] = shared_ptr<DFHalfComplex>(new DFHalfComplex(df, rocoeff[0], iocoeff[0], false, l_coord, aa));
  large_half[1] = shared_ptr<DFHalfComplex>(new DFHalfComplex(df, rocoeff[1], iocoeff[1], false, l_coord, bb));

  for (int i = 0; i != 2; ++i) {
    array<shared_ptr<Matrix>, 2> large_data;
    array<shared_ptr<Matrix>, 2> large_data2;
    large_half[i] = shared_ptr<DFHalfComplex>(new DFHalfComplex(df, rocoeff[i], iocoeff[i], false, l_coord, make_pair(i,i)));
    large_data = large_half[i]->compute_large_Jop(df, trocoeff[i], tiocoeff[i]);
    for(int j = i; j != 2; ++j) {
      large_data2 = large_half[i]->form_2index_large_large(large_half[j]);
      //TODO figure out what to do with this damn data
    }
  }

#if 0
    *rdata_[i] += *df->compute_Jop(rhalfbj[i], trocoeff[i]);
    *idata_[i] -= *df->compute_Jop(ihalfbj[i], tiocoeff[i]);
#endif

  // Small Half Transforms; swapped only needs xy, xz, and yz
#if 0

  array<shared_ptr<DFHalfComplex>, 18> small_half;
  const int a_basis = 2;
  const int b_basis = 3;
  pair<const int, const int> xx = make_pair(0,0);
  pair<const int, const int> xy = make_pair(0,1);
  pair<const int, const int> xz = make_pair(0,2);
  pair<const int, const int> yx = make_pair(1,0);
  pair<const int, const int> yy = make_pair(1,1);
  pair<const int, const int> yz = make_pair(1,2);
  pair<const int, const int> zx = make_pair(2,0);
  pair<const int, const int> zy = make_pair(2,1);
  pair<const int, const int> zz = make_pair(2,2);

  // XX
  small_half[0] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[0], rocoeff[a_basis], iocoeff[a_basis], false, xx, aa));
  small_half[1] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[0], rocoeff[b_basis], iocoeff[b_basis], false, xx, bb));

  // XY
  small_half[2] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[a_basis], iocoeff[a_basis], false, xy, aa));
  small_half[3] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[b_basis], iocoeff[b_basis], false, xy, bb));

  // XZ
  small_half[4] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[a_basis], iocoeff[a_basis], false, xz, ab));
  small_half[5] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[b_basis], iocoeff[b_basis], false, xz, ba));

  // YX
  small_half[6] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[a_basis], iocoeff[a_basis], true, yx, aa));
  small_half[7] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[b_basis], iocoeff[b_basis], true, yx, bb));

  // YY
  small_half[8] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[3], rocoeff[a_basis], iocoeff[a_basis], false, yy, aa));
  small_half[9] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[3], rocoeff[b_basis], iocoeff[b_basis], false, yy, bb));

  // YZ
  small_half[10] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[a_basis], iocoeff[a_basis], false, yz, ab));
  small_half[11] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[b_basis], iocoeff[b_basis], false, yz, ba));

  // ZX
  small_half[12] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[a_basis], iocoeff[a_basis], true, zx, ab));
  small_half[13] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[b_basis], iocoeff[b_basis], true, zx, ba));

  // ZY
  small_half[14] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[a_basis], iocoeff[a_basis], true, zy, ab));
  small_half[15] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[b_basis], iocoeff[b_basis], true, zy, ba));

  // ZZ
  small_half[16] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[5], rocoeff[a_basis], iocoeff[a_basis], false, zz, aa));
  small_half[17] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[5], rocoeff[b_basis], iocoeff[b_basis], false, zz, bb));

  //Multiply all the dfhalfcomplexes together
  for (int i = 0; i != small_half.size(); ++i) {
    for (int j = i; j != small_half.size(); ++j) {
      pair<const int, const int> index = make_pair(small_half[i]->basis().second, small_half[j]->basis().second);
      array<shared_ptr<Matrix>, 2> small_data;
      small_data = small_half[i]->form_2index_small(small_half[j]);
      //TODO
      //insert into Larger Matrix Here?
      //figure out duplicate blocks
    }
  }
  
  pair<const int, const int> eri_coord[6] = {xx, xy, xz, yy, yz, zz};
  pair<const int, const int> eri_basis1[6] = {aa, aa, ab, aa, ab, aa};
  pair<const int, const int> eri_basis2[6] = {bb, bb, ba, bb, ba, bb};

  for (int i = 0; i != dfs.size(); ++i) {
    for (int j = i; j != small_half.size(); ++j) {
      array<shared_ptr<Matrix>, 2> small_Jop_data1;
      array<shared_ptr<Matrix>, 2> small_Jop_data2;
      small_Jop_data1 = small_half[j]->compute_small_Jop(dfs[i], trocoeff, tiocoeff, eri_basis1[i], eri_coord[i]); 
      small_Jop_data2 = small_half[j]->compute_small_Jop(dfs[i], trocoeff, tiocoeff, eri_basis2[i], eri_coord[i]); 
    }
  }

  // Small - Large Transforms
  for (int i = 0; i != 2; ++i) {
    for (int j = 0; j != small_half.size(); ++j) {
      array<shared_ptr<Matrix>, 2> large_small_data;
      large_small_data = large_half[i]->form_2index_small(small_half[j]);
      //TODO figure out what to do with this damn data
    }
  }

#endif
}
