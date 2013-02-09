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

  {
    // get individual df dist objects for each block and add df to dfs
    vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
    dfs.push_back(geom_->df());
    bool gaunt = false;

    list<shared_ptr<DFData>> dfdists = make_dfdists(dfs, gaunt);
    list<shared_ptr<DFHalfComplex>> half_complex = make_half_complex(dfdists, rocoeff, iocoeff);

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
          complex<double> fac = conj((*j)->fac() / (*i)->fac());
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


#if 0
  vector<shared_ptr<const DFDist>> dfsl = geom_->dfsl()->split_blocks();
  list<shared_ptr<DFData>> mixed_dfdists = make_mixed_array(dfsl);
  list<shared_ptr<DFHalfComplex>> mixed_complexes = mixed_complex(mixed_dfdists, rocoeff, iocoeff);
  for (auto& i : mixed_dfdists) {
    for (auto& j : mixed_complexes) {
      const int k = j->coeff_matrix();
      add_mixed_Jop_block(j, i, trocoeff[k], tiocoeff[k]);
    }
  }
#endif
#if 0

  // before computing K operators, we factorize mixed_complexes 
  for (auto i = mixed_complexes.begin(); i != mixed_complexes.end(); ++i) {
    for (auto j = i; j != mixed_complexes.end(); ) {
      if (i != j && (*i)->matches((*j))) {
        complex<double> fac = (*i)->factor(*j);
        (*i)->zaxpy(fac, (*j)); 
        j = mixed_complexes.erase(j);
      } else {
        ++j;
      } 
    }
  }
#endif

  // will use the zgemm3m-like algorithm
#if 0
  for (auto& i : mixed_complexes)
    i->set_sum_diff();
#endif

#if 0
  // computing K operators
  for (auto i = mixed_complexes.begin(); i != mixed_complexes.end(); ++i) {
    for (auto j = mixed_complexes.begin(); j != mixed_complexes.end(); ++j) {
      add_mixed_Exop_block(*i, *j, scale_exchange); 
    }
  }

#endif

}


void DFock::add_Jop_block(list<shared_ptr<DFHalfComplex>> dfc, shared_ptr<const DFData> dfdata, list<shared_ptr<const ZMatrix>> cd) { 

  const int n = geom_->nbasis();
  const tuple<int, int, int, int> index = dfdata->compute_index_Jop();

  shared_ptr<ZMatrix> sum = cd.front()->clone();

  auto cditer = cd.begin();
  for (auto& i : dfc) { 
    complex<double> coeff1 = conj(i->fac()) * dfdata->fac();
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
    const double dfac = (real(dfdata->fac()) == 0 ? -1.0 : 1.0);
    shared_ptr<ZMatrix> tjop(new ZMatrix(*dat->transpose() * dfac));
    add_block(n * get<1>(index), n * get<0>(index), n, n, tjop);
    add_block(n * get<3>(index), n * get<2>(index), n, n, tjop);
  }
}


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
//*a *= dfc1->coeff1() * conj(dfc1->coeff2()) * conj(dfc2->coeff1()) * dfc2->coeff2();
  *a *= conj(dfc1->fac()) * dfc2->fac();

  int index0, index1;
  tie(index0, index1) = dfc1->compute_index_Exop(dfc2);

  add_block(n*index0, n*index1, n, n, a);

  if (dfc1 != dfc2) {
    add_block(n*index1, n*index0, n, n, a->transpose_conjg());
  }

}


list<shared_ptr<DFData>> DFock::make_dfdists(vector<shared_ptr<const DFDist>> dfs, bool gaunt) {
  list<shared_ptr<DFData>> dfdists;
  auto k = dfs.begin();
  if (!gaunt) { // Regular DHF
    for (int i = 0; i != 3; ++i) {
      for (int j = i; j != 3; ++j) {
        dfdists.push_back(shared_ptr<DFData>(new DFData(*k++, make_pair(i,j), Comp::L)));
      }
    }
    // large-large
    dfdists.push_back(shared_ptr<DFData>(new DFData(*k++, make_pair(3,3), Comp::L)));
    assert(k == dfs.end());
  } else { // Gaunt Term
    assert(false);
    for (int i = 0; i != 3; ++i) {
      for (int alpha = 0; alpha != 3; ++alpha) {
        dfdists.push_back(shared_ptr<DFData>(new DFData(*k++, make_pair(i,Comp::L), alpha)));
      }
    }
    assert(k == dfs.end());
  }

  return dfdists;
}


list<shared_ptr<DFHalfComplex>> DFock::make_half_complex(list<shared_ptr<DFData>> dfdists, array<shared_ptr<const Matrix>,4> rocoeff, 
                                                     array<shared_ptr<const Matrix>,4> iocoeff) {
  list<shared_ptr<DFHalfComplex>> half_complex; 
  for (auto& i : dfdists) {
    half_complex.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(i, rocoeff, iocoeff)));
    half_complex.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(i->opp(), rocoeff, iocoeff)));
    if (i->cross()) {
      half_complex.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(i->swap(), rocoeff, iocoeff)));
      half_complex.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(i->opp_and_swap(), rocoeff, iocoeff)));
    }
  }
  return half_complex;
  
}
