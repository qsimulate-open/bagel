//
// BAGEL - Parallel electron correlation program.
// Filename: dfock.cc
// Copyright (C) 2013 Matthew Kelley
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
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
#include <src/rel/cdmatrix_drv.h>

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

  driver(rocoeff, iocoeff, trocoeff, tiocoeff, false, false, scale_exchange);
  if (gaunt_) {
    driver(rocoeff, iocoeff, trocoeff, tiocoeff, gaunt_, breit_, scale_exchange);
 // driver(rocoeff, iocoeff, trocoeff, tiocoeff, gaunt_, false, scale_exchange);
  }
}


void DFock::add_Jop_block(shared_ptr<const DFData> dfdata, list<shared_ptr<const CDMatrix>> cd, const double scale, bool gaunt, bool breit) {

  const int n = geom_->nbasis();
  vector<shared_ptr<ZMatrix>> dat = dfdata->compute_Jop(cd);

  //add it twice, once to first basis combo, then once to opposite basis combo
  int j = 0;
  for (auto& i : dfdata->basis()) {
    add_block(i->fac()*scale, n * i->basis(0), n * i->basis(1), n, n, dat[j++]);
  }

  //if basis1 != basis2, get transpose to fill in opposite corner
  if (dfdata->cross()) {
    shared_ptr<const DFData> swap = dfdata->swap();
    int j = 0;
    for (auto& i : swap->basis()) {
      add_block(i->fac()*scale, n * i->basis(0), n * i->basis(1), n, n, dat[j++]->transpose());
    }
  }
}


void DFock::add_Exop_block(shared_ptr<DFHalfComplex> dfc1, shared_ptr<DFHalfComplex> dfc2, const double scale, const bool diag, const bool notranspose) {

  // minus from -1 in the definition of exchange
  const int n = geom_->nbasis();

  shared_ptr<Matrix> r, i;
  if (!dfc1->sum()) {
    cout << "** warning : using 4 multiplication" << endl;
    // plus
    r   =  dfc1->get_real()->form_2index(dfc2->get_real(), 1.0);
    // plus = minus * minux. (one from i*i, the other from conjugate)
    *r += *dfc1->get_imag()->form_2index(dfc2->get_imag(), 1.0);
    // minus (from conjugate)
    i   =  dfc1->get_real()->form_2index(dfc2->get_imag(), -1.0);
    // plus
    *i += *dfc1->get_imag()->form_2index(dfc2->get_real(), 1.0);
  } else {
    // the same as above
    shared_ptr<Matrix> ss = dfc1->sum()->form_2index(dfc2->sum(), 0.5);
    shared_ptr<Matrix> dd = dfc1->diff()->form_2index(dfc2->diff(), 0.5);
    r = shared_ptr<Matrix>(new Matrix(*ss + *dd));
    i = shared_ptr<Matrix>(new Matrix(*ss - *dd + *dfc1->get_real()->form_2index(dfc2->get_imag(), -2.0)));
  }

  const bool diagonal = diag || dfc1 == dfc2;

  shared_ptr<ZMatrix> a(new ZMatrix(*r, *i));
  for (auto& i1 : dfc1->basis()) {
    for (auto& i2 : dfc2->basis()) {
      shared_ptr<ZMatrix> out(new ZMatrix(*a * (conj(i1->fac())*i2->fac())));

      const int index0 = i1->basis(1);
      const int index1 = i2->basis(1);

      add_block(-scale, n*index0, n*index1, n, n, out);

      if ((!diagonal || i1 != i2) && !notranspose) {
        add_block(-scale, n*index1, n*index0, n, n, out->transpose_conjg());
      }
    }
  }

}


list<shared_ptr<DFData>> DFock::make_dfdists(vector<shared_ptr<const DFDist>> dfs, bool mixed) {
  const vector<int> xyz = { Comp::X, Comp::Y, Comp::Z };

  list<shared_ptr<DFData>> dfdists;
  if (!mixed) { // Regular DHF
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
      dfdists.push_back(shared_ptr<DFData>(new DFData(*k++, make_pair(i,Comp::L), xyz)));
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

void DFock::driver(array<shared_ptr<const Matrix>, 4> rocoeff, array<shared_ptr<const Matrix>, 4> iocoeff,
                              array<shared_ptr<const Matrix>, 4> trocoeff, array<shared_ptr<const Matrix>, 4>tiocoeff, bool gaunt, bool breit,
                              const double scale_exchange)  {

  Timer timer(0);

  const double gscale = gaunt ? (breit ? -0.5 : -1.0) : 1.0;

  if (breit && !gaunt)
    throw logic_error("What are you smoking, son?! Don't call breit without gaunt!");

  vector<shared_ptr<const DFDist>> dfs;
  if (!gaunt) {
    // get individual df dist objects for each block and add df to dfs
    dfs = geom_->dfs()->split_blocks();
    dfs.push_back(geom_->df());
  } else if (gaunt) {
    dfs = geom_->dfsl()->split_blocks();
  }

  list<shared_ptr<DFData>> dfdists = make_dfdists(dfs, gaunt);
  list<shared_ptr<DFHalfComplex>> half_complex = make_half_complex(dfdists, rocoeff, iocoeff);

  const string printtag = !gaunt ? "Coulomb" : "Gaunt";
  timer.tick_print(printtag + ": half trans");

  // split
  list<shared_ptr<DFHalfComplex>> half_complex_exch, half_complex_exch2;
  for (auto& i : half_complex) {
    list<shared_ptr<DFHalfComplex>> tmp = i->split(!breit);
    half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
  }

  // before computing K operators, we factorize half_complex
  for (auto i = half_complex_exch.begin(); i != half_complex_exch.end(); ++i) {
    for (auto j = i; j != half_complex_exch.end(); ) {
        if (i != j && (*i)->matches((*j)) && (*i)->alpha_matches((*j))) {
        complex<double> fac = conj((*j)->fac() / (*i)->fac());
        (*i)->zaxpy(fac, (*j));
        j = half_complex_exch.erase(j);
      } else {
        ++j;
      }
    }
  }
  //assert(half_complex_exch.size() == 8);

  if (breit) {
    // first make a copy of half_complex_exch
    for (auto& i : half_complex_exch)
      half_complex_exch2.push_back(i->copy());

    shared_ptr<Breit> breit_matrix(new Breit(geom_));
    list<shared_ptr<Breit2Index>> breit_2index;
    for (int i = 0; i != breit_matrix->nblocks(); ++i) {
      breit_2index.push_back(shared_ptr<Breit2Index>(new Breit2Index(breit_matrix->index(i), breit_matrix->data(i), geom_->df()->data2())));

      // if breit index is xy, xz, yz, get yx, zx, zy (which is the exact same with reversed index)
      if (breit_matrix->cross(i))
        breit_2index.push_back(breit_2index.back()->cross());
    }

    for (auto& i : half_complex_exch) {
      for (auto& j : breit_2index) {
        if (i->alpha_matches(j))
          half_complex_exch2.push_back(i->multiply_breit2index(j));
      }
    }
    timer.tick_print("Breit: 2-index mulitplied");

    for (auto i = half_complex_exch2.begin(); i != half_complex_exch2.end(); ++i) {
      for (auto j = i; j != half_complex_exch2.end(); ) {
          if (i != j && (*i)->matches((*j)) && (*i)->alpha_matches((*j))) {
          complex<double> fac = conj((*j)->fac() / (*i)->fac());
          (*i)->zaxpy(fac, (*j));
          j = half_complex_exch2.erase(j);
        } else {
          ++j;
        }
      }
    }
  } else {
    half_complex_exch2 = half_complex_exch;
  }

  // will use the zgemm3m-like algorithm
  for (auto& i : half_complex_exch)
    i->set_sum_diff();
  if (half_complex_exch != half_complex_exch2) {
    for (auto& i : half_complex_exch2)
      i->set_sum_diff();
  }

  // computing K operators
  int icnt = 0;
  for (auto i = half_complex_exch.begin(); i != half_complex_exch.end(); ++i, ++icnt) {
    int jcnt = 0;
    for (auto j = half_complex_exch2.begin(); j != half_complex_exch2.end(); ++j, ++jcnt) {
      if ((*i)->alpha_matches((*j)) && icnt <= jcnt) {
        add_Exop_block(*i, *j, gscale*scale_exchange, icnt == jcnt);
      }
    }
  }

  timer.tick_print(printtag + ": K operator");

  list<shared_ptr<const CDMatrix>> cd;
  // compute J operators
  for (auto& j : half_complex_exch2) {
    for (auto& i : j->basis()) {
      cd.push_back(shared_ptr<CDMatrix>(new CDMatrix_drv(j, i, trocoeff, tiocoeff, geom_->df()->data2())));
    }
  }

  for (auto& i : dfdists) {
    add_Jop_block(i, cd, gscale, gaunt, breit);
  }

  timer.tick_print(printtag + ": J operator");
}
