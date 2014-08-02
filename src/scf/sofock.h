//
// BAGEL - Parallel electron correlation program.
// Filename: sofock.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __BAGEL_SRC_SCF_SOFOCK_H
#define __BAGEL_SRC_SCF_SOFOCK_H

#include <src/df/df.h>
#include <src/integral/libint/libint.h>
#include <src/math/zmatrix.h>

namespace bagel {

class SOFock : public ZMatrix {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    const std::shared_ptr<const ZMatrix> previous_;
    const std::shared_ptr<const ZMatrix> socoeff_;
    void form_sofock();

  public:
    SOFock(const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const ZMatrix> previous,
           const std::shared_ptr<const ZMatrix> socoeff) :
      ZMatrix(geom->nbasis() * 2, geom->nbasis() * 2), geom_(geom), previous_(previous), socoeff_(socoeff)
    {
      form_sofock();
    }

};

void SOFock::form_sofock() {
  *this += *previous_;

  int const nocc = socoeff_->mdim();
  int const nbasis = socoeff_->ndim() / 2;

  std::shared_ptr<const ZMatrix> coeffa = socoeff_->get_submatrix(0, 0, nbasis, nocc);
  std::shared_ptr<const ZMatrix> coeffb = socoeff_->get_submatrix(nbasis, 0, nbasis, nocc);

  auto rcoeffa = std::make_shared<const Matrix>(*coeffa->get_real_part());
  auto icoeffa = std::make_shared<const Matrix>(*coeffa->get_imag_part());
  auto rcoeffb = std::make_shared<const Matrix>(*coeffb->get_real_part());
  auto icoeffb = std::make_shared<const Matrix>(*coeffb->get_imag_part());

  std::shared_ptr<const DFDist> df = geom_->df();

  std::shared_ptr<DFHalfDist> rhalf_a = df->compute_half_transform(rcoeffa)->apply_J();
  std::shared_ptr<DFHalfDist> ihalf_a = df->compute_half_transform(icoeffa)->apply_J();

  std::shared_ptr<DFHalfDist> rhalf_b = df->compute_half_transform(rcoeffb)->apply_J();
  std::shared_ptr<DFHalfDist> ihalf_b = df->compute_half_transform(icoeffb)->apply_J();

  std::shared_ptr<Matrix> J = std::make_shared<Matrix>(nbasis, nbasis);
  *J += *df->compute_Jop(rhalf_a, rcoeffa->transpose(), true);
  *J += *df->compute_Jop(ihalf_a, icoeffa->transpose(), true);
  *J += *df->compute_Jop(rhalf_b, rcoeffb->transpose(), true);
  *J += *df->compute_Jop(ihalf_b, icoeffb->transpose(), true);

  std::shared_ptr<Matrix> rfockaa = std::make_shared<Matrix>(nbasis, nbasis);
  std::shared_ptr<Matrix> ifockaa = std::make_shared<Matrix>(nbasis, nbasis);
  *rfockaa += *rhalf_a->form_2index(rhalf_a, -1.0);
  *rfockaa += *ihalf_a->form_2index(ihalf_a, -1.0);
  *ifockaa += *ihalf_a->form_2index(rhalf_a, -1.0);
  *ifockaa += *rhalf_a->form_2index(ihalf_a,  1.0);

  std::shared_ptr<Matrix> rfockab = std::make_shared<Matrix>(nbasis, nbasis);
  std::shared_ptr<Matrix> ifockab = std::make_shared<Matrix>(nbasis, nbasis);
  *rfockab += *rhalf_a->form_2index(rhalf_b, -1.0);
  *rfockab += *ihalf_a->form_2index(ihalf_b, -1.0);
  *ifockab += *ihalf_a->form_2index(rhalf_b, -1.0);
  *ifockab += *rhalf_a->form_2index(ihalf_b,  1.0);

  const std::complex<double> real(1.0, 0.0);
  const std::complex<double> imag(0.0, 1.0);

  add_real_block(real, 0, 0, nbasis, nbasis, rfockaa);
  add_real_block(imag, 0, 0, nbasis, nbasis, ifockaa);
  add_real_block(real, 0, 0, nbasis, nbasis, J);
  add_real_block(real, nbasis, nbasis, nbasis, nbasis, rfockaa);
  add_real_block(-imag, nbasis, nbasis, nbasis, nbasis, ifockaa);
  add_real_block(real, nbasis, nbasis, nbasis, nbasis, J);

  add_real_block(real, 0, nbasis, nbasis, nbasis, rfockab);
  add_real_block(imag, 0, nbasis, nbasis, nbasis, ifockab);
  add_real_block(real, nbasis, 0, nbasis, nbasis, rfockab->transpose());
  add_real_block(-imag, nbasis, 0, nbasis, nbasis, ifockab->transpose());
}

}

#endif
