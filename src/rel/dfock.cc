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

  // get individual df dist objects for each block and add df to dfs
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  vector<shared_ptr<const DFDist>> dfsl = geom_->dfsl()->split_blocks();

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

  list<shared_ptr<DFHalfComplex>> half_complex;
  list<shared_ptr<DFData>> dfdists;
  tie(half_complex, dfdists) = make_arrays(rocoeff, iocoeff, dfs);

  // compute J operators
  list<shared_ptr<const ZMatrix>> cd;
  for (auto& j : half_complex) {
    const int k =  j->coeff_matrix();
    cd.push_back(shared_ptr<ZMatrix>(new ZMatrix(
     *j->get_real()->compute_cd(trocoeff[k], geom_->df()->data2(), true)+*j->get_imag()->compute_cd(tiocoeff[k], geom_->df()->data2(), true),
     *j->get_real()->compute_cd(tiocoeff[k], geom_->df()->data2(), true)-*j->get_imag()->compute_cd(trocoeff[k], geom_->df()->data2(), true))));
  }
  for (auto& i : dfdists) {
    add_Jop_block(half_complex, i, cd); 
  }

  // before computing K operators, we factorize half_complex 
  for (auto i = half_complex.begin(); i != half_complex.end(); ++i) {
    for (auto j = i; j != half_complex.end(); ) {
      if (i != j && (*i)->matches((*j))) {
        complex<double> fac = (*i)->factor(*j);
        (*i)->zaxpy(fac, (*j)); 
        j = half_complex.erase(j);
      } else {
        ++j;
      } 
    }
  }
  assert(half_complex.size() == 8);

  // will use the zgemm3m-like algorithm
  for (auto& i : half_complex)
    i->set_sum_diff();

  // computing K operators
  for (auto i = half_complex.begin(); i != half_complex.end(); ++i) {
    for (auto j = i; j != half_complex.end(); ++j) {
      add_Exop_block(*i, *j, scale_exchange); 
    }
  }

}


void DFock::add_Jop_block(list<shared_ptr<DFHalfComplex>> dfc, shared_ptr<const DFData> dfdata, list<shared_ptr<const ZMatrix>> cd) { 

  const int n = geom_->nbasis();
  const tuple<int, int, int, int> index = dfdata->compute_index_Jop();

  shared_ptr<ZMatrix> sum = cd.front()->clone();

  auto cditer = cd.begin();
  for (auto& i : dfc) { 
    complex<double> coeff1 = i->coeff1() * conj(i->coeff2()) * conj(dfdata->coeff1()) * dfdata->coeff2();
    sum->zaxpy(coeff1, *cditer++);
  }

  shared_ptr<Matrix> rdat = dfdata->df()->compute_Jop_from_cd(sum->get_real_part());
  shared_ptr<Matrix> idat = dfdata->df()->compute_Jop_from_cd(sum->get_imag_part());
  shared_ptr<const ZMatrix> dat(new ZMatrix(*rdat, *idat));

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

#if 0
void DFock::add_mixed_Jop_block(list<shared_ptr<DFHalfComplex>> dfsl, shared_ptr<const DFData> dfdata, list<shared_ptr<const ZMatrix>> cd) { 

  const int n = geom_->nbasis();
  const tuple<int, int, int, int> index = dfdata->compute_index_mixed_Jop();

  shared_ptr<ZMatrix> sum = cd.front()->clone();

  auto cditer = cd.begin();
  for (auto& i : dfc) { 
    complex<double> coeff1 = i->compute_coeff(dfdata);
    sum->zaxpy(coeff1, *cditer++);
  }

  shared_ptr<Matrix> rdat = dfdata->df()->compute_Jop_from_cd(sum->get_real_part());
  shared_ptr<Matrix> idat = dfdata->df()->compute_Jop_from_cd(sum->get_imag_part());
  shared_ptr<const ZMatrix> dat(new ZMatrix(*rdat, *idat));

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
#endif


void DFock::add_Exop_block(shared_ptr<DFHalfComplex> dfc1, shared_ptr<DFHalfComplex> dfc2, const double scale) {

  // minus from -1 in the definition of exchange
  const double scale_exch = - scale;
  const int n = geom_->nbasis();

  shared_ptr<Matrix> r, i;
  if (!dfc1->sum()) {
    cout << "** warning : using 4 multiplication" << endl;
    // plus
    r   =  dfc1->get_real()->form_2index(dfc2->get_real(), scale_exch); 
    // plus = minus * minux. (one from i*i, the other from conjugate)
    *r += *dfc1->get_imag()->form_2index(dfc2->get_imag(), scale_exch);
    // minus (from conjugate)
    i   =  dfc1->get_real()->form_2index(dfc2->get_imag(), -scale_exch);
    // plus
    *i += *dfc1->get_imag()->form_2index(dfc2->get_real(), scale_exch);
  } else {
    // the same as above
    shared_ptr<Matrix> ss = dfc1->sum()->form_2index(dfc2->sum(), -0.5);
    shared_ptr<Matrix> dd = dfc1->diff()->form_2index(dfc2->diff(), -0.5);
    r = shared_ptr<Matrix>(new Matrix(*ss + *dd));
    i = shared_ptr<Matrix>(new Matrix(*ss - *dd + *dfc1->get_real()->form_2index(dfc2->get_imag(), 2.0)));
  }

  shared_ptr<ZMatrix> a(new ZMatrix(*r, *i));
  *a *= dfc1->coeff1() * conj(dfc1->coeff2()) * conj(dfc2->coeff1()) * dfc2->coeff2();

  int index0, index1;
  tie(index0, index1) = dfc1->compute_index_Exop(dfc2);

  add_block(n*index0, n*index1, n, n, a);

  if (dfc1 != dfc2) {
    add_block(n*index1, n*index0, n, n, a->transpose_conjg());
  }

}


tuple<list<shared_ptr<DFHalfComplex>>, list<shared_ptr<DFData>>>
DFock::make_arrays(array<shared_ptr<const Matrix>,4> rocoeff, array<shared_ptr<const Matrix>,4> iocoeff, vector<shared_ptr<const DFDist>> dfs) {

  list<shared_ptr<DFData>> dfdists;
  auto k = dfs.begin();
  for (int i = 0; i != 3; ++i) {
    for (int j = i; j != 3; ++j) {
      dfdists.push_back(shared_ptr<DFData>(new DFData(*k++, make_pair(i,j))));
    }
  }
  // large-large
  dfdists.push_back(shared_ptr<DFData>(new DFData(*k++, make_pair(3,3))));
  assert(k == dfs.end());

  list<shared_ptr<DFHalfComplex>> half_complex; 
  for (auto& i : dfdists) {
    const int coeff_index = i->coeff_index();
    half_complex.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(i,              rocoeff[coeff_index],   iocoeff[coeff_index])));
    half_complex.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(i->opp(),       rocoeff[coeff_index+1], iocoeff[coeff_index+1])));
    if (i->cross()) {
      half_complex.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(i->swap(),         rocoeff[coeff_index],   iocoeff[coeff_index])));
      half_complex.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(i->opp_and_swap(), rocoeff[coeff_index+1], iocoeff[coeff_index+1])));
    }
  }

  return make_tuple(half_complex, dfdists);
}

list<shared_ptr<DFData>> DFock::make_mixed(vector<shared_ptr<const DFDist>> dfsl) {
  list<shared_ptr<DFData>> mixed_dfdists;
  auto k = dfsl.begin();
  for (int i = 0; i != 3; ++i) {
    mixed_dfdists.push_back(shared_ptr<DFData>(new DFData(*k, make_pair(i,DFData::Comp::L))));
    mixed_dfdists.push_back(shared_ptr<DFData>(new DFData(*k++, make_pair(DFData::Comp::L,i))));
  }
  assert(k == dfsl.end());
  return mixed_dfdists;
}

list<shared_ptr<DFHalfComplex>> DFock::mixed_complex(list<shared_ptr<DFData>> mixed_dfdists, array<shared_ptr<const Matrix>,4> rocoeff, 
                                                     array<shared_ptr<const Matrix>,4> iocoeff) {
  list <shared_ptr<DFHalfComplex>> mixed_complex;
  for (auto& i : mixed_dfdists) {
    //TODO figure out index
    const int coeff_index = i->coeff_index();
    mixed_complex.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(i,        rocoeff[coeff_index], iocoeff[coeff_index])));
    mixed_complex.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(i->opp(), rocoeff[coeff_index+1], iocoeff[coeff_index+1])));
  }

  return mixed_complex;
}

