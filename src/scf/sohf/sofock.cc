//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: sofock.cc
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


#include <src/util/math/matrix.h>
#include <src/scf/sohf/sofock.h>

using namespace bagel;
using namespace std;

SOFock::SOFock(const shared_ptr<const Geometry> geom, const shared_ptr<const ZMatrix> previous, const shared_ptr<const ZMatrix> socoeff) :
  ZMatrix(geom->nbasis() * 2, geom->nbasis() * 2), geom_(geom), previous_(previous), socoeff_(socoeff)
{ form_sofock(); }

void SOFock::form_sofock() {
  *this += *previous_;

  int const nocc = socoeff_->mdim();
  int const nbasis = socoeff_->ndim() / 2;

  shared_ptr<const ZMatrix> coeffa = socoeff_->get_submatrix(0, 0, nbasis, nocc);
  shared_ptr<const ZMatrix> coeffb = socoeff_->get_submatrix(nbasis, 0, nbasis, nocc);

  auto rcoeffa = make_shared<const Matrix>(*coeffa->get_real_part());
  auto icoeffa = make_shared<const Matrix>(*coeffa->get_imag_part());
  auto rcoeffb = make_shared<const Matrix>(*coeffb->get_real_part());
  auto icoeffb = make_shared<const Matrix>(*coeffb->get_imag_part());

  shared_ptr<const DFDist> df = geom_->df();

  shared_ptr<DFHalfDist> rhalf_a = df->compute_half_transform(rcoeffa)->apply_J();
  shared_ptr<DFHalfDist> ihalf_a = df->compute_half_transform(icoeffa)->apply_J();

  shared_ptr<DFHalfDist> rhalf_b = df->compute_half_transform(rcoeffb)->apply_J();
  shared_ptr<DFHalfDist> ihalf_b = df->compute_half_transform(icoeffb)->apply_J();

  shared_ptr<Matrix> J = make_shared<Matrix>(nbasis, nbasis);
  *J += *df->compute_Jop(rhalf_a, rcoeffa->transpose(), true);
  *J += *df->compute_Jop(ihalf_a, icoeffa->transpose(), true);
  *J += *df->compute_Jop(rhalf_b, rcoeffb->transpose(), true);
  *J += *df->compute_Jop(ihalf_b, icoeffb->transpose(), true);

  shared_ptr<Matrix> rfockaa = make_shared<Matrix>(nbasis, nbasis);
  shared_ptr<Matrix> ifockaa = make_shared<Matrix>(nbasis, nbasis);
  *rfockaa += *rhalf_a->form_2index(rhalf_a, -1.0);
  *rfockaa += *ihalf_a->form_2index(ihalf_a, -1.0);
  *ifockaa += *ihalf_a->form_2index(rhalf_a, -1.0);
  *ifockaa += *rhalf_a->form_2index(ihalf_a,  1.0);

  shared_ptr<Matrix> rfockab = make_shared<Matrix>(nbasis, nbasis);
  shared_ptr<Matrix> ifockab = make_shared<Matrix>(nbasis, nbasis);
  *rfockab += *rhalf_a->form_2index(rhalf_b, -1.0);
  *rfockab += *ihalf_a->form_2index(ihalf_b, -1.0);
  *ifockab += *ihalf_a->form_2index(rhalf_b, -1.0);
  *ifockab += *rhalf_a->form_2index(ihalf_b,  1.0);

  const complex<double> real(1.0, 0.0);
  const complex<double> imag(0.0, 1.0);

  add_real_block(real, 0, 0, nbasis, nbasis, *rfockaa);
  add_real_block(imag, 0, 0, nbasis, nbasis, *ifockaa);
  add_real_block(real, 0, 0, nbasis, nbasis, *J);
  add_real_block(real, nbasis, nbasis, nbasis, nbasis, *rfockaa);
  add_real_block(-imag, nbasis, nbasis, nbasis, nbasis, *ifockaa);
  add_real_block(real, nbasis, nbasis, nbasis, nbasis, *J);

  add_real_block(real, 0, nbasis, nbasis, nbasis, *rfockab);
  add_real_block(imag, 0, nbasis, nbasis, nbasis, *ifockab);
  add_real_block(real, nbasis, 0, nbasis, nbasis, *rfockab->transpose());
  add_real_block(-imag, nbasis, 0, nbasis, nbasis, *ifockab->transpose());
}
