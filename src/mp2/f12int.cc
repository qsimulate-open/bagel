//
// Newint - Parallel electron correlation program.
// Filename: f12int.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/mp2/f12int.h>
#include <iostream>

using namespace std;


F12Int::F12Int(const multimap<string, string> id, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re, const double gam)
 : idata_(id), geom_(geom), ref_(re), gamma_(gam) {

  // somewhat naive implementations based on 4-index MO integrals *incore*

  // coefficient sets
  const size_t nocc = geom->nocc();
  const size_t nbasis = geom->nbasis();
  const size_t nvirt = nbasis - nocc;
  const double* const oc = ref_->coeff()->data();
  const double* const vc = oc + nocc*nbasis; 
 
  shared_ptr<F12Mat> ymat;
  {
  // Yukawa integral can be thrown right away
  shared_ptr<YukawaFit> yukawa(new YukawaFit(geom->nbasis(), geom->naux(), gamma_,
                               geom->atoms(), geom->offsets(), geom->aux_atoms(), geom->aux_offsets(),
                               0.0, geom->df()));

  const shared_ptr<const DF_Full> yoo = yukawa->compute_half_transform(oc, nocc)->compute_second_transform(oc, nocc);
  shared_ptr<F12Mat> ym(new F12Mat(nocc, yoo->form_4index(yoo))); ymat = ym;
  }

  shared_ptr<SlaterFit> slater(new SlaterFit(geom->nbasis(), geom->naux(), gamma_,
                               geom->atoms(), geom->offsets(), geom->aux_atoms(), geom->aux_offsets(),
                               0.0, geom->df()));
  
};
