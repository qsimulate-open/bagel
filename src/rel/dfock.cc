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

  // get individual df dist objects for each block and add df to dfs
  vector<shared_ptr<const DFDist>> dfs = dfs_total->split_blocks();
  dfs.push_back(df);

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
  vector<shared_ptr<DFData>> dfdists;
  make_arrays(rocoeff, iocoeff, dfs, half_complex, dfdists);

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

void DFock::add_Jop_block(shared_ptr<DFHalfComplex> dfc, shared_ptr<const DFData> dfdata, shared_ptr<const Matrix> trocoeff, 
                 shared_ptr<const Matrix> tiocoeff) {

  const int n = geom_->nbasis();
  shared_ptr<const DFDist> df = dfdata->df();
  
  complex<double> coeff1 = dfc->compute_coeff(dfdata->basis(), dfdata->coord());
  complex<double> im(0.0, 1.0);
  double coeff2 = (coeff1.real() == 0.0 ? -1.0 :1.0);
  const tuple<int, int, int, int> index = dfc->compute_index_Jop(dfdata->basis(), dfdata->coord());

  shared_ptr<Matrix> real, imag;
  real   =  df->compute_Jop(dfc->get_real(), trocoeff, true);
  *real += *df->compute_Jop(dfc->get_imag(), tiocoeff, true);
  imag   =  df->compute_Jop(dfc->get_real(), tiocoeff, true);
  *imag -= *df->compute_Jop(dfc->get_imag(), trocoeff, true);
  *imag *= -coeff2;
  shared_ptr<ZMatrix> dat(new ZMatrix(n, n));
  dat->add_real_block(coeff1,    0, 0, n, n, real); 
  dat->add_real_block(coeff1*im, 0, 0, n, n, imag); 

  //add it twice, once to first basis combo, then once to opposite basis combo
  add_block(n * get<0>(index), n * get<1>(index), n, n, dat);
  add_block(n * get<2>(index), n * get<3>(index), n, n, dat);

  //if basis1 != basis2, get transpose to fill in opposite corner
  if (dfdata->cross()) {
    shared_ptr<ZMatrix> tjop(new ZMatrix(*dat->transpose() * dfdata->cross_coeff()));
    add_block(n * get<1>(index), n * get<0>(index), n, n, tjop);
    add_block(n * get<3>(index), n * get<2>(index), n, n, tjop);
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

void DFock::make_arrays(array<shared_ptr<const Matrix>, 4> rocoeff, array<shared_ptr<const Matrix>, 4> iocoeff,
                                 vector<shared_ptr<const DFDist> > dfs, array<shared_ptr<DFHalfComplex>, 20>& half_complex, vector<shared_ptr<DFData>> dfdists) {

  for (int i = 0, k = 0; i != 4; ++i) {
    for (int j = i; j != 4; ++j) {
      if (j == 3 && i != 3) break;
      pair<const int, const int> coord = make_pair(i,j);
      shared_ptr<DFData> dfcoord(new DFData(dfs[k], coord));
      dfdists.push_back(dfcoord);
      k++;
    }
  }

  for (int i = 0, j = 0; i != dfdists.size(); ++i, ++j) {
    int coeff_index = dfdists[i]->coeff_index();
    half_complex[j] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfdists[i], rocoeff[coeff_index], iocoeff[coeff_index]));
    half_complex[++j] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfdists[i]->opp_basis(), rocoeff[coeff_index+1], iocoeff[coeff_index+1]));
    if (dfdists[i]->cross()) {
      half_complex[++j] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfdists[i]->swap_df(), rocoeff[coeff_index], iocoeff[coeff_index]));
      half_complex[++j] = shared_ptr<DFHalfComplex>(new DFHalfComplex(dfdists[i]->opp_and_swap(), rocoeff[coeff_index+1], iocoeff[coeff_index+1]));
    }
  }

}

