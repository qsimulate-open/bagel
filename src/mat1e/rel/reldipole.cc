//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reldipole.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/mat1e/rel/reldipole.h>
#include <src/mat1e/dipolematrix.h>
#include <src/mat1e/rel/small1e.h>
#include <src/integral/os/mmbatch.h>

using namespace std;
using namespace bagel;

array<shared_ptr<const ZMatrix>,3> RelDipole::compute_matrices() const {
  // first calculate AO integrals.
  shared_ptr<const Matrix1eArray<3>>  large = make_shared<const DipoleMatrix>(geom_);
  shared_ptr<const Matrix1eArray<12>> small = make_shared<const Small1e<DipoleBatch, shared_ptr<const Molecule>>>(geom_, geom_);

  const int n = geom_->nbasis();

  array<shared_ptr<const ZMatrix>,3> out;

  // loop over cartesian components
  for (int x = 0; x != 3; ++x) {
    auto zdip = make_shared<ZMatrix>(2*n, 2*n);
    zdip->copy_real_block(1.0, 0, 0, n, n, (*large)[x]);
    zdip->copy_real_block(1.0, n, n, n, n, (*large)[x]);

    // inefficient.. (but probably it does not matter)
    const complex<double> w(0.25/(c__*c__));
    const complex<double> wi(0,w.real());
    auto smalldip = zdip->clone();
    smalldip->add_real_block(w,   0, 0, n, n,  (*small)[0+x*4]);
    smalldip->add_real_block(w,   n, n, n, n,  (*small)[0+x*4]);
    smalldip->add_real_block(wi,  0, 0, n, n,  (*small)[1+x*4]);
    smalldip->add_real_block(-wi, n, n, n, n,  (*small)[1+x*4]);
    smalldip->add_real_block(wi,  0, n, n, n,  (*small)[2+x*4]);
    smalldip->add_real_block(wi,  n, 0, n, n,  (*small)[2+x*4]);
    smalldip->add_real_block(w,   0, n, n, n,  (*small)[3+x*4]);
    smalldip->add_real_block(-w,  n, 0, n, n,  (*small)[3+x*4]);

    auto work = make_shared<ZMatrix>(n*4, n*4);
    work->copy_block(0, 0, 2*n, 2*n, zdip);
    work->copy_block(2*n, 2*n, 2*n, 2*n, smalldip);
    out[x] = work;
  }
  return out;
}


array<double,3> RelDipole::compute() const {
  if (density_ == nullptr)
    throw logic_error("RelDipole::compute was called without density matrix");

  array<double,3> out;
  array<shared_ptr<const ZMatrix>,3> matrices = compute_matrices();
  for (int x = 0; x != 3; ++x) {
    out[x] = density_->dot_product(matrices[x]).real();
    assert(fabs(density_->dot_product(matrices[x]).imag()) < 1.0e-6);
  }

  cout << "    * Permanent dipole moment:" << (jobname_.empty() ? "" : " " + jobname_) << endl;
  cout << "           (" << setw(12) << setprecision(6) << out[0] << ", " << setw(12) << out[1] << ", " << setw(12) << out[2] << ") a.u." << endl;

  return out;
}
