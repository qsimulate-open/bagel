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

  array<shared_ptr<DFHalfComplex>, 20> half_complex;
  array<shared_ptr<DFData>, 7> dfdists;
  make_arrays(rocoeff, iocoeff, df, dfs, half_complex, dfdists);

  for (auto& i : dfdists) {
    for (auto& j : half_complex) {
      const int k =  j->coeff_matrix();
      add_Jop_block(j, i, trocoeff[k], tiocoeff[k]); 
    }
  }

  for (auto& i : half_complex)
    i->set_sum_diff();

  for (auto i = half_complex.begin(); i != half_complex.end(); ++i) {
    for (auto j = i; j != half_complex.end(); ++j) {
      add_Exop_block(*i, *j); 
    }
  }

}

void DFock::make_arrays(array<shared_ptr<const Matrix>, 4> rocoeff, array<shared_ptr<const Matrix>, 4> iocoeff, shared_ptr<const DFDist> df, 
                                 vector<shared_ptr<DFDist> > dfs, array<shared_ptr<DFHalfComplex>, 20>& half_complex, array<shared_ptr<DFData>, 7>& dfdists) {

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

  //large
  dfdists[0] = shared_ptr<DFData>(new DFData(df, l_coord, aa));
  //small
  dfdists[1] = shared_ptr<DFData>(new DFData(dfs[0], xx, aa));
  dfdists[2] = shared_ptr<DFData>(new DFData(dfs[1], xy, aa));
  dfdists[3] = shared_ptr<DFData>(new DFData(dfs[2], xz, ab));
  dfdists[4] = shared_ptr<DFData>(new DFData(dfs[3], yy, aa));
  dfdists[5] = shared_ptr<DFData>(new DFData(dfs[4], yz, ab));
  dfdists[6] = shared_ptr<DFData>(new DFData(dfs[5], zz, aa));

  //large
  half_complex[0] = shared_ptr<DFHalfComplex>(new DFHalfComplex(df, rocoeff[0], iocoeff[0], false, l_coord, aa));
  half_complex[1] = shared_ptr<DFHalfComplex>(new DFHalfComplex(df, rocoeff[1], iocoeff[1], false, l_coord, bb));

  // XX
  half_complex[2] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[0], rocoeff[2], iocoeff[2], false, xx, aa));
  half_complex[3] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[0], rocoeff[3], iocoeff[3], false, xx, bb));

  // XY
  half_complex[4] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[2], iocoeff[2], false, xy, aa));
  half_complex[5] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[3], iocoeff[3], false, xy, bb));

  // XZ
  half_complex[6] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[2], iocoeff[2], false, xz, ab));
  half_complex[7] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[3], iocoeff[3], false, xz, ba));

  // YX
  half_complex[8] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[2], iocoeff[2], true, yx, aa));
  half_complex[9] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[1], rocoeff[3], iocoeff[3], true, yx, bb));

  // YY
  half_complex[10] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[3], rocoeff[2], iocoeff[2], false, yy, aa));
  half_complex[11] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[3], rocoeff[3], iocoeff[3], false, yy, bb));

  // YZ
  half_complex[12] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[2], iocoeff[2], false, yz, ab));
  half_complex[13] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[3], iocoeff[3], false, yz, ba));

  // ZX
  half_complex[14] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[2], iocoeff[2], true, zx, ab));
  half_complex[15] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[2], rocoeff[3], iocoeff[3], true, zx, ba));

  // ZY
  half_complex[16] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[2], iocoeff[2], true, zy, ab));
  half_complex[17] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[4], rocoeff[3], iocoeff[3], true, zy, ba));

  // ZZ
  half_complex[18] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[5], rocoeff[2], iocoeff[2], false, zz, aa));
  half_complex[19] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfs[5], rocoeff[3], iocoeff[3], false, zz, bb));
}

void DFock::add_Jop_block(shared_ptr<DFHalfComplex> dfc, shared_ptr<const DFData> dfdata, shared_ptr<const Matrix> trocoeff, 
                 shared_ptr<const Matrix> tiocoeff) {

  const int n = geom_->nbasis();
  shared_ptr<const DFDist> df = dfdata->df();
  
  complex<double> coeff1 = dfc->compute_coeff(dfdata->basis(), dfdata->coord());
  double coeff2 = (coeff1.real() == 0.0 ? -1.0 :1.0);
  const tuple<int, int, int, int> index = dfc->compute_index_Jop(dfdata->basis(), dfdata->coord());

  array<shared_ptr<Matrix>, 4> jop;
  jop[0] = df->compute_Jop(dfc->get_real(), trocoeff, true);
  jop[1] = df->compute_Jop(dfc->get_imag(), tiocoeff, true);
  jop[2] = df->compute_Jop(dfc->get_real(), tiocoeff, true);
  jop[3] = df->compute_Jop(dfc->get_imag(), trocoeff, true);
  *jop[2] *= -coeff2;
  *jop[3] *= coeff2;

  for (int i = 0; i != 4; ++i) {
   //add it twice, once to first basis combo, then once to opposite basis combo
    add_real_block(coeff1, n * get<0>(index), n * get<1>(index), n, n, jop[i]);
    add_real_block(coeff1, n * get<2>(index), n * get<3>(index), n, n, jop[i]);

    //if basis1 != basis2, get transpose to fill in opposite corner
    if (dfdata->cross()) {
      shared_ptr<Matrix> tjop = jop[i]->transpose();
      add_real_block(coeff1, n * get<1>(index), n * get<0>(index), n, n, (*tjop * dfdata->cross_coeff()).data());
      add_real_block(coeff1, n * get<3>(index), n * get<2>(index), n, n, (*tjop * dfdata->cross_coeff()).data());
    }
  }
}

void DFock::add_Exop_block(shared_ptr<DFHalfComplex> dfc1, shared_ptr<DFHalfComplex> dfc2) {

  const int n = geom_->nbasis();
  complex<double> coeff1 = dfc1->compute_coeff(dfc2->basis(), dfc2->coord());
  const tuple<int, int, int, int> index = dfc1->compute_index_Exop(dfc2->basis(), dfc2->coord());

  shared_ptr<Matrix> r, i;
  if (!dfc1->sum()) {
    cout << "** warning : using 4 multiplication" << endl;
    r   =  dfc1->get_real()->form_2index(dfc2->get_real(), -1.0); 
    *r += *dfc1->get_imag()->form_2index(dfc2->get_imag(), -1.0);
    i   =  dfc1->get_real()->form_2index(dfc2->get_imag(),  1.0);
    *i += *dfc1->get_imag()->form_2index(dfc2->get_real(), -1.0);
  } else {
    shared_ptr<Matrix> ss = dfc1->sum()->form_2index(dfc2->sum(), -0.5);
    shared_ptr<Matrix> dd = dfc1->diff()->form_2index(dfc2->diff(), -0.5);
    r = shared_ptr<Matrix>(new Matrix(*ss + *dd));
    i = shared_ptr<Matrix>(new Matrix(*ss - *dd + *dfc1->get_real()->form_2index(dfc2->get_imag(), 2.0)));
  }

  shared_ptr<ZMatrix> a(new ZMatrix(r->ndim(), r->mdim()));
  a->add_real_block(coeff1, 0, 0, n, n, r);
  complex<double> im(0.0, 1.0);
  double coeff2 = (coeff1.real() == 0.0 ? -1.0 :1.0);
  a->add_real_block(coeff1*im*coeff2, 0, 0, n, n, i);

  add_block(n * get<0>(index), n * get<1>(index), n, n, a);

  if (dfc1 != dfc2) {
    add_block(n * get<1>(index), n * get<0>(index), n, n, a->transpose_conjg());
  }

}

