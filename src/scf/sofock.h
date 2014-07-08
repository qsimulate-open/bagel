//
// BAGEL - Parallel electron correlation program.
// Filename: sofock.h
// Copyright (C) 2012 Toru Shiozaki
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
//#include <src/integral/rys/eribatch.h>

namespace bagel {

class SOFock : public Matrix {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    const std::shared_ptr<const Matrix> previous_;
    const std::shared_ptr<const Matrix> socoeff_;
    void form_sofock();

  public:
    SOFock(const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Matrix> previous,
           const std::shared_ptr<const Matrix> socoeff) : Matrix(geom->nbasis() * 2, geom->nbasis() * 2), geom_(geom), previous_(previous),
           socoeff_(socoeff) {
      form_sofock();
    }

};

void SOFock::form_sofock() {
  *this += *previous_;

  int const nocc = socoeff_->mdim();
  int const nbasis = socoeff_->ndim() / 2;

  std::shared_ptr<const Matrix> coeffa = socoeff_->get_submatrix(0, 0, nbasis, nocc);
  std::shared_ptr<const Matrix> coeffb = socoeff_->get_submatrix(nbasis, 0, nbasis, nocc);

  auto coeffta = std::make_shared<const Matrix>(*coeffa->transpose());
  auto coefftb = std::make_shared<const Matrix>(*coeffb->transpose());

  std::shared_ptr<const DFDist> df = geom_->df();

  std::shared_ptr<DFHalfDist> half_a = df->compute_half_transform(coeffa)->apply_J();
  std::shared_ptr<DFHalfDist> half_b = df->compute_half_transform(coeffb)->apply_J();

  std::shared_ptr<Matrix> fockaa = std::make_shared<Matrix>(nbasis, nbasis);
  std::shared_ptr<Matrix> fockbb = std::make_shared<Matrix>(nbasis, nbasis);
  std::shared_ptr<Matrix> fockab = std::make_shared<Matrix>(nbasis, nbasis);
  std::shared_ptr<Matrix> fockba = std::make_shared<Matrix>(nbasis, nbasis);

  *fockaa += *half_a->form_2index(half_a, -1.0);
  *fockaa += *df->compute_Jop(half_a, coeffta, true);
  *fockbb += *df->compute_Jop(half_a, coeffta, true);

  *fockbb += *half_b->form_2index(half_b, -1.0);
  *fockbb += *df->compute_Jop(half_b, coefftb, true);
  *fockaa += *df->compute_Jop(half_b, coefftb, true);

  *fockab += *half_a->form_2index(half_b, -1.0);
  *fockba += *fockab->transpose();
   
  add_block(1.0, 0, 0, nbasis, nbasis, fockaa);
  add_block(1.0, nbasis, nbasis, nbasis, nbasis, fockbb);
  add_block(1.0, 0, nbasis, nbasis, nbasis, fockab);
  add_block(1.0, nbasis, 0, nbasis, nbasis, fockba);
}

}

#endif
