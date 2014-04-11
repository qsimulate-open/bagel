//
// BAGEL - Parallel electron correlation program.
// Filename: debug_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
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

#include <src/london/debug_london.h>
#include <src/london/get_complexeri.h>

using namespace bagel;
using namespace std;

BOOST_CLASS_EXPORT_IMPLEMENT(Debug_London)

static double get_average (shared_ptr<ZMatrix> S) {
  complex<double>* data = S->data();
  double out = 0.0;
  const int size = S->ndim() * S->mdim();
  for (int i=0; i!=size; i++) out += abs(data[i]);
  out /= size;
  return out;
}


static double get_average (shared_ptr<Matrix> S) {
  double* data = S->data();
  double out = 0.0;
  const int size = S->ndim() * S->mdim();
  for (int i=0; i!=size; i++) out += abs(data[i]);
  out /= size;
  return out;
}


static double get_max (shared_ptr<ZMatrix> S) {
  complex<double>* data = S->data();
  double out = 0.0;
  const int size = S->ndim() * S->mdim();
  for (int i=0; i!=size; i++) {
    double mag = abs(data[i]);
    if (mag > out) out = mag;
  }
  return out;
}


static double get_max (shared_ptr<Matrix> S) {
  double* data = S->data();
  double out = 0.0;
  const int size = S->ndim() * S->mdim();
  for (int i=0; i!=size; i++) {
    double mag = abs(data[i]);
    if (mag > out) out = mag;
  }
  return out;
}


static pair<int,int> get_position (shared_ptr<ZMatrix> S) {
  complex<double>* data = S->data();
  double out = 0.0;
  int pos = 0;
  const int size = S->ndim() * S->mdim();
  for (int i=0; i!=size; i++) {
    double mag = abs(data[i]);
    if (mag > out) {
      out = mag;
      pos = i;
    }
  }
  pair <int, int> position (pos/S->ndim(), pos%S->ndim());
  return position;
}


static pair<int,int> get_position (shared_ptr<Matrix> S) {
  double* data = S->data();
  double out = 0.0;
  int pos = 0;
  const int size = S->ndim() * S->mdim();
  for (int i=0; i!=size; i++) {
    double mag = abs(data[i]);
    if (mag > out) {
      out = mag;
      pos = i;
    }
  }
  pair <int, int> position (pos/S->ndim(), pos%S->ndim());
  return position;
}


Debug_London::Debug_London(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry_London> cgeom, const std::shared_ptr<const Reference> re)
        : Method (idata, cgeom, re) { }

void Debug_London::compute() {
  cout << "  Running debug code for London orbital density fitting." << endl;

  {
    shared_ptr<const DFDist_London> df = cgeom_->df();
    shared_ptr<const ZMatrix> ij = df->data2();

    const int asize = df->block(0)->asize();
    const int b1size = df->block(0)->b1size();
    const int b2size = df->block(0)->b2size();
    shared_ptr<const ZMatrix> jcd   = df->block(0)->get_block(0,asize,0,b1size,0,b2size);
    //shared_ptr<const ZMatrix> abit  = df->block(0)->get_block_conj(0,asize,0,b1size,0,b2size);
    shared_ptr<const ZMatrix> jcds = make_shared<const ZMatrix>(*ij * *jcd);
    //shared_ptr<const ZMatrix> abits = make_shared<const ZMatrix>(*(ij->transpose()) * *abit);
    //shared_ptr<const ZMatrix> abis  = abits->transpose();
    shared_ptr<const ZMatrix> abis  = jcds->transpose();

//    ij->print("half-inverted 2-index matrix", 40);
//    jcd->transpose()->print("3-index matrix (ab|i)", 40);
//    abit->transpose()->print("3-index matrix (ab|i*)", 40);

    const shared_ptr<ZMatrix> ERI = get_ERI(cgeom_);
    const shared_ptr<ZMatrix> DFERI = make_shared<ZMatrix>(*abis * *jcds);
    const shared_ptr<ZMatrix> DFERROR = make_shared<ZMatrix>(*ERI - *DFERI);

    DFERI->print("4-index matrix by density fitting", 30);
    ERI->print("Analytical London ERI Matrix!", 30);
    DFERROR->print("Errors of density fitting", 30);
    cout << setprecision(9) << endl;
    cout << " Average error = " << get_average (DFERROR) << endl;
    cout << "     Max error = " << get_max (DFERROR) << endl;
    pair <int, int> position =  get_position (DFERROR);
    cout << "   At position = " << position.first+1 << ", " << position.second+1 << endl << endl;

    cout << "reference value = " << ERI->element(position.first, position.second) << endl;
    cout << "   fitted value = " << DFERI->element(position.first, position.second) << endl;
    cout << "error of result = " << DFERROR->element(position.first, position.second) << endl;
  }

  {
    shared_ptr<const DFDist> df = cgeom_->form_fit<DFDist_ints<ERIBatch>>(cgeom_->overlap_thresh(), true);

    const int asize = df->block(0)->asize();
    const int b1size = df->block(0)->b1size();
    const int b2size = df->block(0)->b2size();
    shared_ptr<const Matrix> jcd = df->get_block(0,asize,0,b1size,0,b2size);
    shared_ptr<const Matrix> ij = df->data2();

    shared_ptr<const Matrix> jcds = make_shared<const Matrix>(*ij * *jcd);
    shared_ptr<const Matrix> abi = jcd->transpose();
    shared_ptr<const Matrix> abis = jcds->transpose();

//    ij->print("half-inverted 2-index matrix", 20);
//    abi->print("3-index matrix (ab|i)", 20);

    const shared_ptr<Matrix> ERI = get_ERI_original(cgeom_);
    const shared_ptr<Matrix> DFERI = make_shared<Matrix>(*abis * *jcds);
    const shared_ptr<Matrix> DFERROR = make_shared<Matrix>(*ERI - *DFERI);

    DFERI->print("4-index matrix by density fitting", 30);
    ERI->print("Analytical London ERI Matrix!", 30);
    DFERROR->print("Errors of density fitting", 30);
    cout << setprecision(9) << endl;
    cout << " Average error = " << get_average (DFERROR) << endl;
    cout << "     Max error = " << get_max (DFERROR) << endl;
    pair <int, int> position =  get_position (DFERROR);
    cout << "   At position = " << position.first+1 << ", " << position.second+1 << endl << endl;

    cout << "reference value = " << ERI->element(position.first, position.second) << endl;
    cout << "   fitted value = " << DFERI->element(position.first, position.second) << endl;
    cout << "error of result = " << DFERROR->element(position.first, position.second) << endl;
  }
}


shared_ptr<const Reference> Debug_London::conv_to_ref() const {
//  auto out = make_shared<Reference>(geom_, coeff(), nocc(), 0, coeff_->mdim()-nocc(), energy());
//  out->set_eig(eig_);
//  return out;
  return nullptr;
}
