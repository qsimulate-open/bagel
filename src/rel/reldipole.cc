//
// BAGEL - Parallel electron correlation program.
// Filename: reldipole.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/rel/reldipole.h>
#include <src/molecule/dipolematrix.h>
#include <src/integral/os/mmbatch.h>
#include <src/rel/small1e.h>

using namespace std;
using namespace bagel;

array<double,3> RelDipole::compute() {
  // first calculate AO integrals.
  shared_ptr<const Matrix1eArray<3>> large = make_shared<const DipoleMatrix>(geom_);
  shared_ptr<const Matrix1eArray<12>> small = make_shared<const Small1e<DipoleBatch>>(geom_);

  const int n = geom_->nbasis();

  array<double,3> out;

  // loop over cartesian components
  for (int x = 0; x != 3; ++x) {
    auto zdip = make_shared<ZMatrix>(2*n, 2*n);
    zdip->copy_real_block(1.0, 0, 0, n, n, (*large)[x]);
    zdip->copy_real_block(1.0, n, n, n, n, (*large)[x]);

    // inefficient.. (but probably it does not matter)
    const complex<double> w(0.25/(c__*c__));
    const complex<double> wi(0,w.real());
    array<shared_ptr<ZMatrix>,4> zsmalldip;
    for (auto& i : zsmalldip)
      i = zdip->clone();
    zsmalldip[0]->copy_real_block(w,   0, 0, n, n,  (*small)[0+x*4]);
    zsmalldip[0]->copy_real_block(w,   n, n, n, n,  (*small)[0+x*4]);
    zsmalldip[1]->copy_real_block(wi,  0, 0, n, n,  (*small)[1+x*4]);
    zsmalldip[1]->copy_real_block(-wi, n, n, n, n,  (*small)[1+x*4]);
    zsmalldip[2]->copy_real_block(wi,  0, n, n, n,  (*small)[2+x*4]);
    zsmalldip[2]->copy_real_block(wi,  n, 0, n, n,  (*small)[2+x*4]);
    zsmalldip[3]->copy_real_block(w,   0, n, n, n,  (*small)[3+x*4]);
    zsmalldip[3]->copy_real_block(-w,  n, 0, n, n,  (*small)[3+x*4]);

    auto smalldip = make_shared<ZMatrix>(*zsmalldip[0] + *zsmalldip[1] + *zsmalldip[2] + *zsmalldip[3]);

    auto data = make_shared<ZMatrix>(n*4, n*4);
    data->copy_block(0, 0, 2*n, 2*n, zdip);
    data->copy_block(2*n, 2*n, 2*n, 2*n, smalldip);
    out[x] = density_->dot_product(data).real();
    assert(fabs(density_->dot_product(data).imag()) < 1.0e-6);
  }

  cout << "    * Permanent dipole moment:" << (jobname_.empty() ? "" : " " + jobname_) << endl;
  cout << "           (" << setw(12) << setprecision(6) << out[0] << ", " << setw(12) << out[1] << ", " << setw(12) << out[2] << ") a.u." << endl;

  return out;
}
