//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: sohcore.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/mat1e/sohcore.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(SOHcore)

SOHcore::SOHcore(shared_ptr<const Molecule> geom, shared_ptr<const Hcore> h)
            : ZMatrix(2*geom->nbasis(), 2*geom->nbasis()), hcore_(h) {
  form_sohcore();
}

void SOHcore::form_sohcore() {

  const int nbasis = ndim()/2;
  const complex<double> real(1.0, 0.0);
  const complex<double> imag(0.0, 1.0);

  hcore_->hso()->allreduce();
  hcore_->hso()->fill_upper_negative();

  add_real_block(real, 0, 0, nbasis, nbasis, *hcore_);
  add_real_block(imag, 0, 0, nbasis, nbasis, *hcore_->hso()->soiaa());

  add_real_block(real, nbasis, nbasis, nbasis, nbasis, *hcore_);
  add_real_block(-imag, nbasis, nbasis, nbasis, nbasis, *hcore_->hso()->soiaa());

  add_real_block(real, 0, nbasis, nbasis, nbasis, *hcore_->hso()->sorab());
  add_real_block(imag, 0, nbasis, nbasis, nbasis, *hcore_->hso()->soiab());

  add_real_block(-real, nbasis, 0, nbasis, nbasis, *hcore_->hso()->sorab());
  add_real_block( imag, nbasis, 0, nbasis, nbasis, *hcore_->hso()->soiab());

}
