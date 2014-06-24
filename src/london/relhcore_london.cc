//
// BAGEL - Parallel electron correlation program.
// Filename: relhcore_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/util/constants.h>
#include <src/london/relhcore_london.h>
#include <src/integral/compos/complexoverlapbatch.h>

using namespace std;
using namespace bagel;

void RelHcore_London::compute_() {
  const int n = geom_->nbasis();

  auto nai     = make_shared<ZMatrix>(2*n, 2*n);
  nai->copy_block(0, 0, n, n, nai_);
  nai->copy_block(n, n, n, n, nai_);

  auto kinetic = make_shared<ZMatrix>(2*n, 2*n);
  kinetic->copy_block(0, 0, n, n, kinetic_);
  kinetic->copy_block(n, n, n, n, kinetic_);

  const complex<double> w(0.25/(c__*c__));
  const complex<double> wi(0.0, w.real());
  auto zsnai = make_shared<ZMatrix>(2*n, 2*n);
  zsnai->add_block(  w, 0, 0, n, n, (*smallnai_)[0]);
  zsnai->add_block(  w, n, n, n, n, (*smallnai_)[0]);
  zsnai->add_block( wi, 0, 0, n, n, (*smallnai_)[1]);
  zsnai->add_block(-wi, n, n, n, n, (*smallnai_)[1]);
  zsnai->add_block( wi, 0, n, n, n, (*smallnai_)[2]);
  zsnai->add_block( wi, n, 0, n, n, (*smallnai_)[2]);
  zsnai->add_block(  w, 0, n, n, n, (*smallnai_)[3]);
  zsnai->add_block( -w, n, 0, n, n, (*smallnai_)[3]);

  // TODO If careful, we should be able to get smalloverlap[0] from ComplexKineticBatch
  //      and smalloverlap[1-3] from ComplexOverlapBatch & magnetic field
  const complex<double> rh (-0.5);
  const complex<double> ih (0.0, rh.real());
  auto zeeman = make_shared<ZMatrix>(2*n, 2*n);
  zeeman->zero();

  //(*(*smalloverlap_)[0]*0.5 - *kinetic_).print("Difference between Kinetic integral and SmallOverlap mimicry", 100);
  //(*(*smalloverlap_)[1]*0.5 - *overlap_*ih*geom_->magnetic_field(2)).print("Difference between Bz-scaled overlap integral and SmallOverlap mimicry", 100);
  //(*(*smalloverlap_)[2]*0.5 - *overlap_*ih*geom_->magnetic_field(0)).print("Difference between Bx-scaled overlap integral and SmallOverlap mimicry", 100);
  //(*(*smalloverlap_)[3]*0.5 - *overlap_*ih*geom_->magnetic_field(1)).print("Difference between By-scaled overlap integral and SmallOverlap mimicry", 100);
  //kinetic_->print("Kinetic energy Matrix", 100);
  //overlap_->print("Overlap Matrix", 100);
  //(*(*smalloverlap_)[0]*0.5).print("Small integral overlap (Zeeman?) Matrix - part 0 (xx + yy + zz)", 100);
  //(*(*smalloverlap_)[1]*0.5).print("Small integral overlap (Zeeman?) Matrix - part 1 (xy - yx)", 100);
  //(*(*smalloverlap_)[2]*0.5).print("Small integral overlap (Zeeman?) Matrix - part 2 (yz - zy)", 100);
  //(*(*smalloverlap_)[3]*0.5).print("Small integral overlap (Zeeman?) Matrix - part 3 (zx - xz)", 100);

  //zeeman->add_block( rh, 0, 0, n, n, (*smalloverlap_)[0]);
  //zeeman->add_block( rh, n, n, n, n, (*smalloverlap_)[0]);
  {
    /*
    zeeman->add_block(  ih, 0, 0, n, n, (*smalloverlap_)[1]);
    zeeman->add_block( -ih, n, n, n, n, (*smalloverlap_)[1]);
    zeeman->add_block( -ih, 0, n, n, n, (*smalloverlap_)[2]);
    zeeman->add_block( -ih, n, 0, n, n, (*smalloverlap_)[2]);
    zeeman->add_block(  rh, 0, n, n, n, (*smalloverlap_)[3]);
    zeeman->add_block( -rh, n, 0, n, n, (*smalloverlap_)[3]);
    */
  }
  {
    ///*
    const complex<double> r2 (0.5);
    const complex<double> i2 (0.0, r2.real());
    zeeman->add_block( -r2*geom_->magnetic_field(2), 0, 0, n, n, overlap_);
    zeeman->add_block(  r2*geom_->magnetic_field(2), n, n, n, n, overlap_);
    zeeman->add_block( -r2*geom_->magnetic_field(0), 0, n, n, n, overlap_);
    zeeman->add_block( -r2*geom_->magnetic_field(0), n, 0, n, n, overlap_);
    zeeman->add_block(  i2*geom_->magnetic_field(1), 0, n, n, n, overlap_);
    zeeman->add_block( -i2*geom_->magnetic_field(1), n, 0, n, n, overlap_);
    //*/
  }

  // RKB hcore: T is off diagonal block matrices, V is first main diagonal, and 1/4m^2c^2W-T is second main diagonal
  zero();
  copy_block(0,   0, 2*n, 2*n, nai);
  copy_block(0, 2*n, 2*n, 2*n, kinetic);
  copy_block(2*n, 0, 2*n, 2*n, kinetic);
  copy_block(2*n, 2*n, 2*n, 2*n, zsnai);
  add_block(-1.0, 2*n, 2*n, 2*n, 2*n, kinetic);
  add_block( -1.0, 0, 2*n, 2*n, 2*n, zeeman);
  add_block( -1.0, 2*n, 0, 2*n, 2*n, zeeman);
  add_block(  1.0, 2*n, 2*n, 2*n, 2*n, zeeman);

}

