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
#include <src/math/matrix.h>
#include <src/rel/cdmatrix.h>

using namespace std;
using namespace bagel;

// TODO batch size should be automatically determined by the memory size etc.
const static int batchsize = 250;

void DFock::two_electron_part(const shared_ptr<const ZMatrix> coeff, const bool rhf, const double scale_exchange) {

  if (!rhf) throw logic_error("DFock::two_electron_part() is not implemented for non RHF cases");
  assert(geom_->nbasis()*4 == coeff->ndim());

  auto ocoeffall = make_shared<ZMatrix>(*coeff);
  const int nocc = coeff->mdim();
  const int nbatch = (nocc-1) / batchsize+1;
  StaticDist dist(nocc, nbatch);
  vector<pair<size_t, size_t>> table = dist.atable();

  for (auto& itable : table) {
    // Separate Coefficients into real and imaginary
    array<shared_ptr<const Matrix>, 4> rocoeff;
    array<shared_ptr<const Matrix>, 4> iocoeff;
    array<shared_ptr<const Matrix>, 4> trocoeff;
    array<shared_ptr<const Matrix>, 4> tiocoeff;

    for (int i = 0; i != 4; ++i) {
      shared_ptr<const ZMatrix> ocoeff = coeff->get_submatrix(i*geom_->nbasis(), itable.first, geom_->nbasis(), itable.second);
      rocoeff[i] = ocoeff->get_real_part();
      iocoeff[i] = ocoeff->get_imag_part();
      trocoeff[i] = rocoeff[i]->transpose();
      tiocoeff[i] = iocoeff[i]->transpose();
    }

    driver(rocoeff, iocoeff, trocoeff, tiocoeff, false, false, scale_exchange);
    if (gaunt_) {
      driver(rocoeff, iocoeff, trocoeff, tiocoeff, gaunt_, breit_, scale_exchange);
    }
  }
}


void DFock::add_Jop_block(shared_ptr<const RelDF> dfdata, list<shared_ptr<const CDMatrix>> cd, const double scale) {

  const int n = geom_->nbasis();
  vector<shared_ptr<ZMatrix>> dat = dfdata->compute_Jop(cd);

  //add it twice, once to first basis combo, then once to opposite basis combo
  int j = 0;
  for (auto& i : dfdata->basis()) {
    add_block(i->fac(dfdata->cartesian())*scale, n * i->basis(0), n * i->basis(1), n, n, dat[j++]);
  }

  //if basis1 != basis2, get transpose to fill in opposite corner
  if (dfdata->not_diagonal()) {
    shared_ptr<const RelDF> swap = dfdata->swap();
    int j = 0;
    for (auto& i : swap->basis()) {
      add_block(i->fac(swap->cartesian())*scale, n*i->basis(0), n*i->basis(1), n, n, dat[j++]->transpose_conjg());
      // conjg does not matter because (1) offdiagonal of Coulomb is real; (2) offdiagonal of Gaunt and Breit is zero. 
      // (2) is due to [sigma_w, sigma_w']_+ = \delta_ww' -- tricky! 
    }
  }
}


void DFock::add_Exop_block(shared_ptr<RelDFHalf> dfc1, shared_ptr<RelDFHalf> dfc2, const double scale, const bool diag) {

  // minus from -1 in the definition of exchange
  const int n = geom_->nbasis();

  shared_ptr<Matrix> r, i;
  if (!dfc1->sum()) {
    cout << "** warning : using 4 multiplication" << endl;
    r   =  dfc1->get_real()->form_2index(dfc2->get_real(), 1.0);
    *r += *dfc1->get_imag()->form_2index(dfc2->get_imag(), 1.0);
    i   =  dfc1->get_real()->form_2index(dfc2->get_imag(), 1.0);
    *i += *dfc1->get_imag()->form_2index(dfc2->get_real(),-1.0);
  } else {
    // the same as above
    shared_ptr<Matrix> ss = dfc1->sum()->form_2index(dfc2->sum(), 0.5);
    shared_ptr<Matrix> dd = dfc1->diff()->form_2index(dfc2->diff(), 0.5);
    r = make_shared<Matrix>(*ss + *dd);
    i = make_shared<Matrix>(*dd - *ss + *dfc1->get_real()->form_2index(dfc2->get_imag(), 2.0));
  }

  const bool diagonal = diag || dfc1 == dfc2;

  auto a = make_shared<ZMatrix>(*r, *i);
  for (auto& i1 : dfc1->basis()) {
    for (auto& i2 : dfc2->basis()) {
      auto out = make_shared<ZMatrix>(*a * (conj(i1->fac(dfc1->cartesian()))*i2->fac(dfc2->cartesian())));

      const int index0 = i1->basis(1);
      const int index1 = i2->basis(1);

      add_block(-scale, n*index0, n*index1, n, n, out);

      if (!diagonal || *i1 != *i2) {
        add_block(-scale, n*index1, n*index0, n, n, out->transpose_conjg());
      }
    }
  }

}


list<shared_ptr<RelDF>> DFock::make_dfdists(vector<shared_ptr<const DFDist>> dfs, bool mixed) {
  const vector<int> xyz = { Comp::X, Comp::Y, Comp::Z };

  list<shared_ptr<RelDF>> dfdists;
  if (!mixed) { // Regular DHF
    const vector<int> alphaL = { Comp::L };
    auto k = dfs.begin();
    for (auto& i : xyz) {
      for (auto& j : xyz)
        if (i <= j)
          dfdists.push_back(make_shared<RelDF>(*k++, make_pair(i,j), vector<int>{Comp::L}));
    }
    // large-large
    dfdists.push_back(make_shared<RelDF>(*k++, make_pair(Comp::L,Comp::L), vector<int>{Comp::L}));
    assert(k == dfs.end());

  } else { // Gaunt Term
    auto k = dfs.begin();
    for (auto& i : xyz)
      dfdists.push_back(make_shared<RelDF>(*k++, make_pair(i,Comp::L), xyz));
    assert(k == dfs.end());
  }
  return dfdists;
}


list<shared_ptr<RelDFHalf>> DFock::make_half_complex(list<shared_ptr<RelDF>> dfdists, array<shared_ptr<const Matrix>,4> rocoeff,
                                                     array<shared_ptr<const Matrix>,4> iocoeff) {
  list<shared_ptr<RelDFHalf>> half_complex;
  for (auto& i : dfdists) {
    vector<shared_ptr<RelDFHalf>> dat = i->compute_half_transform(rocoeff, iocoeff);
    half_complex.insert(half_complex.end(), dat.begin(), dat.end());

    if (i->not_diagonal()) {
      vector<shared_ptr<RelDFHalf>> dat = i->swap()->compute_half_transform(rocoeff, iocoeff);
      half_complex.insert(half_complex.end(), dat.begin(), dat.end());
    }
  }
  return half_complex;

}


void DFock::driver(array<shared_ptr<const Matrix>, 4> rocoeff, array<shared_ptr<const Matrix>, 4> iocoeff,
                              array<shared_ptr<const Matrix>, 4> trocoeff, array<shared_ptr<const Matrix>, 4>tiocoeff, bool gaunt, bool breit,
                              const double scale_exchange)  {

  Timer timer(0);

  vector<shared_ptr<const DFDist>> dfs;
  if (!gaunt) {
    // get individual df dist objects for each block and add df to dfs
    dfs = geom_->dfs()->split_blocks();
    dfs.push_back(geom_->df());
  } else if (gaunt) {
    dfs = geom_->dfsl()->split_blocks();
  }

  list<shared_ptr<RelDF>> dfdists = make_dfdists(dfs, gaunt);
  // Note that we are NOT using dagger-ed coefficients! -1 factor for imagnary will be compensated by CDMatrix and Exop
  list<shared_ptr<RelDFHalf>> half_complex = make_half_complex(dfdists, rocoeff, iocoeff);

  const string printtag = !gaunt ? "Coulomb" : "Gaunt";
  timer.tick_print(printtag + ": half trans");

  // split
  list<shared_ptr<RelDFHalf>> half_complex_exch, half_complex_exch2;
  for (auto& i : half_complex) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
  }
  half_complex.clear();

  // before computing K operators, we factorize half_complex
  factorize(half_complex_exch);
  assert(gaunt  || half_complex_exch.size() == 8);
  assert(!gaunt || half_complex_exch.size() == 24);

  if (breit) {
    // first make a copy of half_complex_exch
    for (auto& i : half_complex_exch)
      half_complex_exch2.push_back(i->copy());

    auto breitint = make_shared<BreitInt>(geom_);
    list<shared_ptr<Breit2Index>> breit_2index;
    for (int i = 0; i != breitint->nblocks(); ++i) {
      breit_2index.push_back(make_shared<Breit2Index>(breitint->index(i), breitint->data(i), geom_->df()->data2()));

      // if breit index is xy, xz, yz, get yx, zx, zy (which is the exact same with reversed index)
      if (breitint->not_diagonal(i))
        breit_2index.push_back(breit_2index.back()->cross());
    }

    for (auto& i : half_complex_exch) {
      for (auto& j : breit_2index) {
        if (i->alpha_matches(j)) {
          half_complex_exch2.push_back(i->multiply_breit2index(j));
          factorize(half_complex_exch2);
        }
      }
    }
    timer.tick_print("Breit: 2-index mulitplied");

  } else {
    half_complex_exch2 = half_complex_exch;
  }

  // this is a necessary condition if we use symmetry below (Exop) 
  assert(half_complex_exch.size() == half_complex_exch2.size());

  // will use the zgemm3m-like algorithm
  for (auto& i : half_complex_exch)
    i->set_sum_diff();
  if (half_complex_exch != half_complex_exch2) {
    for (auto& i : half_complex_exch2)
      i->set_sum_diff();
  }

  const double gscale = gaunt ? (breit ? -0.5 : -1.0) : 1.0;

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
      cd.push_back(make_shared<CDMatrix>(j, i, trocoeff, tiocoeff, geom_->df()->data2()));
    }
  }

  for (auto& i : dfdists) {
    add_Jop_block(i, cd, gscale);
  }

  timer.tick_print(printtag + ": J operator");
}
