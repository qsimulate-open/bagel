//
// BAGEL - Parallel electron correlation program.
// Filename: complexparalleldf.cc
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

#include <src/df/complexdf_base.h>

using namespace std;
using namespace bagel;


void ComplexDF_base::assign_complex_blocks(ParallelDF& source) {
  const int n = source.block().size();
  if (n % 2 != 0) throw runtime_error("ComplexDF_base requires an even number of blocks in the reference object");
  for (int i=0; i<n; i+=2) {
    real_block_.push_back(source.block(i));
    imag_block_.push_back(source.block(i+1));
  }
  cserial_ = source.serial();
  cnaux_ = source.naux();
  assert(real_block_.size() + imag_block_.size() == n);
}


shared_ptr<ZVectorB> ComplexDF_base::complex_compute_cd(const shared_ptr<const ZMatrix> den, shared_ptr<const Matrix> dat2, const bool onlyonce) const {
  if (real_block_.size() != 1 || imag_block_.size() !=1 ) throw logic_error("compute_Jop so far assumes block_.size() == 2 for complex integrals");
  if (!dat2) throw logic_error("ComplexDF_base::compute_cd was called without 2-index integrals");

  auto outr = make_shared<VectorB>(cnaux_);
  auto outi = make_shared<VectorB>(cnaux_);
  const shared_ptr<Matrix> dr = den->get_real_part();
  const shared_ptr<Matrix> di = den->get_imag_part();

  // TODO Using 4-multiplication
  // D = (D|rs)*d_rs
  shared_ptr<VectorB> tmpr = real_block_[0]->form_vec(dr);
  *tmpr -= *imag_block_[0]->form_vec(di);
  shared_ptr<VectorB> tmpi = real_block_[0]->form_vec(di);
  *tmpi += *imag_block_[0]->form_vec(dr);

  auto tmp = make_shared<ZVectorB>(*tmpr, *tmpi);
  copy_n(tmpr->data(), real_block_[0]->asize(), outr->data()+real_block_[0]->astart());
  copy_n(tmpi->data(), imag_block_[0]->asize(), outi->data()+imag_block_[0]->astart());

  // All reduce
  if (!cserial_) {
    outr->allreduce();
    outi->allreduce();
  }

  *outr = *dat2 * *outr;
  *outi = *dat2 * *outi;
  if (!onlyonce) {
    *outr = *dat2 * *outr;
    *outi = *dat2 * *outi;
  }
  return make_shared<ZVectorB>(*outr, *outi);
}
