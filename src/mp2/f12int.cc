//
// BAGEL - Parallel electron correlation program.
// Filename: f12int.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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


#include <stddef.h>
#include <src/df/df.h>
#include <src/mp2/f12int.h>
#include <iostream>
#include <iomanip>

// TODO this is WRONG!!!!!!!! TODO TODO
#include <src/integral/rys/eribatch.h>
#define SlaterFit DFDist_ints<ERIBatch>
#define YukawaFit DFDist_ints<ERIBatch>

using namespace std;
using namespace bagel;

#if 0
int perm(const int i, const int j, const int k, const int l){
  int retval;
  int ival, jval, kval, lval;
  retval = 1;
  if(j < i) { ival = j; jval = i; }
  else { ival = i; jval = j; }
  if(l < k) { kval = l; lval = k; }
  else { kval = k; lval = l; }
  if(ival != jval) retval = retval*2;
  if(kval != lval) retval = retval*2;
  if((ival != kval) || (jval != lval))retval = retval*2;
  return retval;
}

void print_local(orz::DTensor inp) {
  for (int i=0; i!=OUTSIZE; ++i) {
  for (int j=0; j!=OUTSIZE; ++j) {
  for (int k=0; k!=OUTSIZE; ++k) {
  for (int l=0; l!=OUTSIZE; ++l) {
    cout << fixed << setw(9) << setprecision(6) << inp(l,j,k,i);
  }} cout << endl; }}
}

void print_local1(orz::DTensor inp) {
  cout << endl;
  for (int i=0; i!= inp.shape(0); ++i) {
  for (int j=0; j!= inp.shape(1); ++j) {
    cout << fixed << setw(11) << setprecision(7) << inp(i,j);
  } cout << endl;}
  cout << endl;
}
#endif

static void ddot_to_amp(const F12Mat& inp, int n, double gam, string name = "anon") {
  cout << endl;
  double sing = 1.0/gam/gam*(10.0/64.0);
  double trip = 1.0/gam/gam*(6.0/64.0);
  double a = 0.0;
  for (int i=0; i!=n; ++i) {
  for (int j=0; j!=n; ++j) {
  for (int k=0; k!=n; ++k) {
  for (int l=0; l!=n; ++l) {
    if (l == j && k == i && l == k) {
      a+= inp.data(l,j,k,i)* (sing+trip);
    } else if (l == j && k == i) {
      a+= inp.data(l,j,k,i) * (2.0*sing-trip);
    } else if (l == i && k == j) {
      a+= inp.data(l,j,k,i) * (2.0*trip-sing);
    }
  } } } }
  cout << setw(8) << name << " " << fixed << setw(20) << setprecision(15) << a << endl;
  cout << endl;
}


F12Int::F12Int(const multimap<string, string> id, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re,
               const double gam, const int ncore)
 : idata_(id), geom_(geom), ref_(re), gamma_(gam) {

  // somewhat naive implementations based on 4-index MO integrals *incore*

  // coefficient sets
  const size_t nocc = geom->nele()/2 - ncore;
  const size_t nbasis = geom->nbasis();
  const size_t nvirt = nbasis - nocc - ncore;
  shared_ptr<const Matrix> oc = ref_->coeff()->slice(ncore, ncore+nocc);

  const shared_ptr<const DFDist> df = geom_->df();
  const shared_ptr<const DFHalfDist> dxo = df->compute_half_transform(oc)->apply_J();
  const shared_ptr<const DFFullDist> doo = dxo->compute_second_transform(oc);

  shared_ptr<F12Mat> ymat;
  {
  // Yukawa integral can be thrown right away
  shared_ptr<DFDist> yukawa = geom->form_fit<YukawaFit>(0.0, false, gamma_);
  const shared_ptr<const DFHalfDist> yxo = yukawa->compute_half_transform(oc);
  const shared_ptr<const DFFullDist> yoo = yxo->compute_second_transform(oc)->apply_J(geom->df());
  ymat = robust_fitting(doo, yoo);
  }

  shared_ptr<SlaterFit> slater = geom->form_fit<SlaterFit>(0.0, false, gamma_);

  // debug area
  ddot_to_amp(*ymat, nocc, gamma_, "Y matrix ");
}


// df and slater are supposed to be J-applied.
shared_ptr<F12Mat> F12Int::robust_fitting(shared_ptr<const DFFullDist> doo, shared_ptr<const DFFullDist> yoo) {
  const int nocc = yoo->nocc1();
  assert(nocc == yoo->nocc2());

  shared_ptr<F12Mat> ym0(new F12Mat(nocc, yoo->form_4index(doo, 1.0)));
  ddot_to_amp(*ym0, nocc, gamma_, "Y matrix orig");

  shared_ptr<const DFFullDist> doo_J = doo->apply_J();
  shared_ptr<const DFFullDist> doo_JS = doo_J->apply_J(yoo->df());
  shared_ptr<F12Mat> ym1(new F12Mat(nocc, doo_JS->form_4index(doo_J, 1.0)));
  ddot_to_amp(*ym1, nocc, gamma_, "Y matrix orig");
  shared_ptr<F12Mat> ym(new F12Mat(*ym0*2.0 - *ym1));
  return ym;
}

