//
// BAGEL - Parallel electron correlation program.
// Filename: zcasbfgs.cc
// Copyright (C) 2014 Jefferson E. Bates
//
// Author: Toru Shiozaki <jefferson.bates@northwestern.edu>
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

#include <src/zcasscf/zcasbfgs.h>
#include <src/zcasscf/zqvec.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;


complex<double> ZCASBFGS::find_level_shift(shared_ptr<const ZRotFile> rotmat) const {
  /* returns the smallest Hessian value that is not below -mc^2 to be used as a level shift
     This particular choice of level shift parameter ensures that the initial diagonal guess has Np negative values
     where Np is the number of positronic orbitals */
  double csq = c__*c__;
  complex<double> l0 = rotmat->data(0);

  for (int j = 1; j != rotmat->size(); ++j) {
    if (l0.real() > rotmat->data(j).real() && csq + rotmat->data(j).real() > 0)
      l0 = rotmat->data(j);
  }
  complex<double> level_shift;
  if (real(l0) < 0) {
    double scale = idata_->get<double>("scalefac", 1.05);
    level_shift = l0 * scale;
    cout << " " << endl;
    cout << setprecision(8) << "Smallest Hessian Element (excluding positrons) = " << l0 << endl;
    cout << setprecision(2) << "Scaling Factor                                 = " << scale << endl;
    cout << setprecision(8) << "Level Shift                                    = " << level_shift << endl << endl;
  } else {
    cout << " level shift is not negative " << endl;
    level_shift = complex<double> (0.0,0.0);
  }

  return level_shift;
}


tuple<shared_ptr<ZRotFile>, vector<double>, shared_ptr<ZRotFile>, shared_ptr<ZRotFile>, bool> ZCASBFGS::optimize_subspace_rotations(vector<double> energy, shared_ptr<const ZRotFile> grad, shared_ptr<const ZRotFile> rot, shared_ptr<SRBFGS<ZRotFile>> srbfgs, shared_ptr<ZMatrix> cold, bool optimize_electrons) {
  // function to optimize only a single subspace of orbitals neglecting the coupling to the other
  const int nvirtnr = nvirt_ - nneg_/2;

  // copy the inputs
  shared_ptr<ZRotFile> newgrad;
  shared_ptr<ZRotFile> newrot = rot->copy();
  if (optimize_electrons) {
    newgrad = copy_electronic_rotations(grad);
  } else {
    newgrad = copy_positronic_rotations(grad);
  }

  const bool tight = idata_->get<bool>("tight", false);
  const int limmem = idata_->get<int>("limited_memory", 0);
  auto reset = srbfgs->check_step(energy, newgrad, newrot, tight, limmem);
  if (reset) {
    cout << " STEP DOES NOT MEET PROPER CRITERIA " << endl;
    // TODO : implement step rejection or another clever alternative
  }

  shared_ptr<ZRotFile> a;
  if (optimize_electrons) {
    a = srbfgs->more_sorensen_extrapolate(newgrad, newrot);
  } else {
    // positronic optimization results in a negative level shift so use Hebden method
    a = srbfgs->extrapolate(newgrad, newrot);
  }

  if (optimize_electrons) {
    kramers_adapt(a, nclosed_, nact_, nvirtnr);
  } else {
    kramers_adapt(a, nclosed_, nact_, nneg_/2);
  }
  cout << setprecision(6) << " Subspace gradient rms  = " << newgrad->rms() << endl;

  return make_tuple(a, energy, newgrad, newrot, reset);
}
