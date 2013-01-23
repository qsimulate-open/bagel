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
  shared_ptr<const DFDist> dfs_total = geom_->form_fit<DFDist_ints<SmallERIBatch> >(1.0e-8, false); // TODO thresh should be controlled from the input deck

  // get individual df dist objects for each block
  vector<shared_ptr<DFDist> > dfs = dfs_total->split_blocks();

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
  array<shared_ptr<DFHalfDist>, 2> rhalfbj;
  array<shared_ptr<DFHalfDist>, 2> ihalfbj;
  array<shared_ptr<DFHalfDist>, 2> rhalf;
  array<shared_ptr<DFHalfDist>, 2> ihalf;

  for (int i = 0; i != 2; ++i) {
    rhalfbj[i] = df->compute_half_transform(rocoeff[i]);
    ihalfbj[i] = df->compute_half_transform(iocoeff[i]);

    //TODO Check Dimensions
    rdata_[i] = shared_ptr<Matrix>(new Matrix(rocoeff[i]->ndim(), rocoeff[i]->mdim()));
    idata_[i] = shared_ptr<Matrix>(new Matrix(iocoeff[i]->ndim(), iocoeff[i]->mdim()));
    rdata_[i]->zero();
    idata_[i]->zero();

    rhalf[i] = rhalfbj[i]->apply_J();
    ihalf[i] = ihalfbj[i]->apply_J();

    *rdata_[i] += *rhalf[i]->form_2index(rhalf[i], -1.0*scale_exchange);
    *idata_[i] += *ihalf[i]->form_2index(ihalf[i], -1.0*scale_exchange);

    // Real and imaginary cross term
    shared_ptr<Matrix> cross = rhalf[i]->form_2index(ihalf[i], -2.0*scale_exchange);
    cross->symmetrize();
    *rdata_[i] += *cross;

    *rdata_[i] += *df->compute_Jop(rhalfbj[i], trocoeff[i]);
    *idata_[i] -= *df->compute_Jop(ihalfbj[i], tiocoeff[i]);
  }

  // Small Half Transforms; swapped only needs xy, xz, and yz

  array<shared_ptr<DFHalfComplex>, 18> small_half;
  const int a_basis = 2;
  const int b_basis = 3;

  // XX
  small_half[0] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[0], rocoeff[a_basis], iocoeff[a_basis], false, make_pair(0,0), make_pair(0,0)));
  small_half[1] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[0], rocoeff[b_basis], iocoeff[b_basis], false, make_pair(0,0), make_pair(1,1)));

  // XY
  small_half[2] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[a_basis], iocoeff[a_basis], false, make_pair(0,1), make_pair(0,0)));
  small_half[3] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[b_basis], iocoeff[b_basis], false, make_pair(0,1), make_pair(1,1)));

  // XZ
  small_half[4] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[a_basis], iocoeff[a_basis], false, make_pair(0,2), make_pair(0,1)));
  small_half[5] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[b_basis], iocoeff[b_basis], false, make_pair(0,2), make_pair(1,0)));

  // YX
  small_half[6] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[a_basis], iocoeff[a_basis], true, make_pair(1,0), make_pair(0,0)));
  small_half[7] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[b_basis], iocoeff[b_basis], true, make_pair(1,0), make_pair(1,1)));

  // YY
  small_half[8] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[3], rocoeff[a_basis], iocoeff[a_basis], false, make_pair(1,1), make_pair(0,0)));
  small_half[9] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[3], rocoeff[b_basis], iocoeff[b_basis], false, make_pair(1,1), make_pair(1,1)));

  // YZ
  small_half[10] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[a_basis], iocoeff[a_basis], false, make_pair(1,2), make_pair(0,1)));
  small_half[11] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[b_basis], iocoeff[b_basis], false, make_pair(1,2), make_pair(1,0)));

  // ZX
  small_half[12] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[a_basis], iocoeff[a_basis], true, make_pair(2,0), make_pair(0,1)));
  small_half[13] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[b_basis], iocoeff[b_basis], true, make_pair(2,0), make_pair(1,0)));

  // ZY
  small_half[14] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[a_basis], iocoeff[a_basis], true, make_pair(2,1), make_pair(0,1)));
  small_half[15] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[b_basis], iocoeff[b_basis], true, make_pair(2,1), make_pair(1,0)));

  // ZZ
  small_half[16] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[5], rocoeff[a_basis], iocoeff[a_basis], false, make_pair(2,2), make_pair(0,0)));
  small_half[17] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[5], rocoeff[b_basis], iocoeff[b_basis], false, make_pair(2,2), make_pair(1,1)));

  //Multiply all the dfhalfcomplexes together
  array<shared_ptr<Matrix>, 648> small_k;
  for (int i = 0; i != small_half.size(); ++i) {
    for (int j = i; j != small_half.size(); ++j) {
      pair<const int, const int> index = make_pair(small_half[i]->basis().second, small_half[j]->basis().second);
      array<shared_ptr<Matrix>, 2> small_data;
      small_data = small_half[i]->form_2index_complex(small_half[j]);
      //insert into Larger Matrix Here?
      //figure out duplicate blocks
    }
  }


}
