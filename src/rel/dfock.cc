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
#include <src/util/timer.h>

using namespace std;
using namespace bagel;

void DFock::two_electron_part(const shared_ptr<const ZMatrix> coeff, const bool rhf, const double scale_exchange) {

  if (!rhf) throw logic_error("DFock::two_electron_part() is not implemented for non RHF cases");
  assert(geom_->nbasis()*4 == coeff->ndim());

  // Separate Coefficients into real and imaginary
  array<shared_ptr<const Matrix>, 4> rocoeff;
  array<shared_ptr<const Matrix>, 4> iocoeff;
  array<shared_ptr<const Matrix>, 4> trocoeff;
  array<shared_ptr<const Matrix>, 4> tiocoeff;

  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ocoeff = coeff->get_submatrix(i*geom_->nbasis(), 0, geom_->nbasis(), coeff->mdim()); 
    rocoeff[i] = ocoeff->get_real_part();
    iocoeff[i] = ocoeff->get_imag_part();
    trocoeff[i] = rocoeff[i]->transpose();
    tiocoeff[i] = iocoeff[i]->transpose();
  }

  Timer timer(0);

  {
    // get individual df dist objects for each block and add df to dfs
    vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
    dfs.push_back(geom_->df());
    bool gaunt = false;

    list<shared_ptr<DFData>> dfdists = make_dfdists(dfs, gaunt);
    list<shared_ptr<DFHalfComplex>> half_complex = make_half_complex(dfdists, rocoeff, iocoeff);

    timer.tick_print("Coulomb: half trans");

    // compute J operators
    list<shared_ptr<const ZMatrix>> cd;
    for (auto& j : half_complex) {
      for (auto& i : j->basis()) {
        const int k =  i->basis(1);
        cd.push_back(shared_ptr<ZMatrix>(new ZMatrix(
         *j->get_real()->compute_cd(trocoeff[k], geom_->df()->data2(), true)+*j->get_imag()->compute_cd(tiocoeff[k], geom_->df()->data2(), true),
         *j->get_real()->compute_cd(tiocoeff[k], geom_->df()->data2(), true)-*j->get_imag()->compute_cd(trocoeff[k], geom_->df()->data2(), true))));
      }
    }
    for (auto& i : dfdists) {
      add_Jop_block(half_complex, i, cd); 
    }

    timer.tick_print("Coulomb: J operator");

    // split
    list<shared_ptr<DFHalfComplex>> half_complex_exch;
    for (auto& i : half_complex) {
      list<shared_ptr<DFHalfComplex>> tmp = i->split();
      half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
    }
    // before computing K operators, we factorize half_complex 
    for (auto i = half_complex_exch.begin(); i != half_complex_exch.end(); ++i) {
      for (auto j = i; j != half_complex_exch.end(); ) {
        if (i != j && (*i)->matches((*j))) {
          complex<double> fac = conj((*j)->fac() / (*i)->fac());
          (*i)->zaxpy(fac, (*j)); 
          j = half_complex_exch.erase(j);
        } else {
          ++j;
        } 
      }
    }
    assert(half_complex_exch.size() == 8);

    // will use the zgemm3m-like algorithm
    for (auto& i : half_complex_exch)
      i->set_sum_diff();

    // computing K operators
    for (auto i = half_complex_exch.begin(); i != half_complex_exch.end(); ++i) {
      for (auto j = i; j != half_complex_exch.end(); ++j) {
        add_Exop_block(*i, *j, scale_exchange); 
      }
    }
    timer.tick_print("Coulomb: K operator");
  }


  // Gaunt term
  if (true) {
    vector<shared_ptr<const DFDist>> dfsl = geom_->dfsl()->split_blocks();
    list<shared_ptr<DFData>> mixed_dfdists = make_dfdists(dfsl, true);
    list<shared_ptr<DFHalfComplex>> mixed_complex = make_half_complex(mixed_dfdists, rocoeff, iocoeff);

    timer.tick_print("Gaunt: half trans");

    // compute J operators
    list<shared_ptr<const ZMatrix>> cd;
    for (auto& j : mixed_complex) {
      for (auto& i : j->basis()) {
        const int k =  i->basis(1);
        cd.push_back(shared_ptr<ZMatrix>(new ZMatrix(
         *j->get_real()->compute_cd(trocoeff[k], geom_->df()->data2(), true)+*j->get_imag()->compute_cd(tiocoeff[k], geom_->df()->data2(), true),
         *j->get_real()->compute_cd(tiocoeff[k], geom_->df()->data2(), true)-*j->get_imag()->compute_cd(trocoeff[k], geom_->df()->data2(), true))));
      }
    }
    for (auto& i : mixed_dfdists) {
      add_Jop_block(mixed_complex, i, cd); 
    }

    timer.tick_print("Gaunt: J operator");

    // split
    list<shared_ptr<DFHalfComplex>> half_complex_exch;
    for (auto& i : mixed_complex) {
      list<shared_ptr<DFHalfComplex>> tmp = i->split();
      half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
    }
    // before computing K operators, we factorize mixed_complex 
    for (auto i = half_complex_exch.begin(); i != half_complex_exch.end(); ++i) {
      for (auto j = i; j != half_complex_exch.end(); ) {
        if (i != j && (*i)->matches((*j))) {
          complex<double> fac = conj((*j)->fac() / (*i)->fac());
          (*i)->zaxpy(fac, (*j)); 
          j = half_complex_exch.erase(j);
        } else {
          ++j;
        } 
      }
    }
    assert(half_complex_exch.size() == 8);

    for (auto& i : half_complex_exch)
      i->set_sum_diff();

    // computing K operators
    for (auto i = half_complex_exch.begin(); i != half_complex_exch.end(); ++i) {
      for (auto j = i; j != half_complex_exch.end(); ++j) {
        add_Exop_block(*i, *j, scale_exchange); 
      }
    }

    timer.tick_print("Gaunt: K operator");
  }

}


void DFock::add_Jop_block(list<shared_ptr<DFHalfComplex>> dfc, shared_ptr<const DFData> dfdata, list<shared_ptr<const ZMatrix>> cd) { 

  const int n = geom_->nbasis();

  shared_ptr<ZMatrix> sum = cd.front()->clone();

  auto cditer = cd.begin();
  for (auto& i : dfc) { 
    for (auto& j : i->basis())
      sum->zaxpy(j->fac(), *cditer++);
  }
  assert(cditer == cd.end());

  shared_ptr<Matrix> rdat = dfdata->df()->compute_Jop_from_cd(sum->get_real_part());
  shared_ptr<Matrix> idat = dfdata->df()->compute_Jop_from_cd(sum->get_imag_part());
  shared_ptr<const ZMatrix> dat(new ZMatrix(*rdat, *idat));

  //add it twice, once to first basis combo, then once to opposite basis combo
  for (auto& i : dfdata->basis())
    add_block(n * i->basis(0), n * i->basis(1), n, n, (*dat*i->fac()).data());

  //if basis1 != basis2, get transpose to fill in opposite corner
  if (dfdata->cross()) {
    shared_ptr<const DFData> swap = dfdata->swap();
    shared_ptr<ZMatrix> tdat = dat->transpose();
    for (auto& i : swap->basis())
      add_block(n * i->basis(0), n * i->basis(1), n, n, (*tdat*i->fac()).data());
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
  for (auto& i1 : dfc1->basis()) {
    for (auto& i2 : dfc2->basis()) {
      shared_ptr<ZMatrix> out(new ZMatrix(*a * (conj(i1->fac())*i2->fac())));

      const int index0 = i1->basis(1);
      const int index1 = i2->basis(1);

      add_block(n*index0, n*index1, n, n, out);

      if (dfc1 != dfc2 || i1 != i2) {
        add_block(n*index1, n*index0, n, n, out->transpose_conjg());
      }
    }
  }

}


list<shared_ptr<DFData>> DFock::make_dfdists(vector<shared_ptr<const DFDist>> dfs, bool gaunt) {
  const vector<int> xyz = { Comp::X, Comp::Y, Comp::Z };

  list<shared_ptr<DFData>> dfdists;
  if (!gaunt) { // Regular DHF
    const vector<int> alphaL = { Comp::L };
    auto k = dfs.begin();
    for (auto& i : xyz) {
      for (auto& j : xyz)
        if (i <= j)
          dfdists.push_back(shared_ptr<DFData>(new DFData(*k++, make_pair(i,j), {Comp::L})));
    }
    // large-large
    dfdists.push_back(shared_ptr<DFData>(new DFData(*k++, make_pair(Comp::L,Comp::L), {Comp::L})));
    assert(k == dfs.end());

  } else { // Gaunt Term
    auto k = dfs.begin();
    for (auto& i : xyz)
      dfdists.push_back(shared_ptr<DFData>(new DFData(*k++, make_pair(i,Comp::L), {Comp::X, Comp::Z})));
    assert(k == dfs.end());
  }
  return dfdists;
}


list<shared_ptr<DFHalfComplex>> DFock::make_half_complex(list<shared_ptr<DFData>> dfdists, array<shared_ptr<const Matrix>,4> rocoeff, 
                                                     array<shared_ptr<const Matrix>,4> iocoeff) {
  list<shared_ptr<DFHalfComplex>> half_complex; 
  for (auto& i : dfdists) {
    vector<shared_ptr<DFHalfComplex>> dat = i->compute_half_transform(rocoeff, iocoeff);
    half_complex.insert(half_complex.end(), dat.begin(), dat.end());

    if (i->cross()) {
      vector<shared_ptr<DFHalfComplex>> dat = i->swap()->compute_half_transform(rocoeff, iocoeff);
      half_complex.insert(half_complex.end(), dat.begin(), dat.end());
    }
  }
  return half_complex;
  
}
