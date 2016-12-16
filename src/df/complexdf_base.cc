//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexdf_base.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

  // D = (D|rs)*d_rs  (using zgemm3m-like algorithm for complex multiplication)
  shared_ptr<VectorB> tmpr, tmpi;
  {
    const shared_ptr<Matrix> dr = den->get_real_part();
    const shared_ptr<Matrix> di = den->get_imag_part();
    auto dri = make_shared<Matrix>(*dr + *di);

    tmpr = real_block_[0]->form_vec(dr);
    shared_ptr<VectorB> tmp0 = imag_block_[0]->form_vec(di);

    shared_ptr<DFBlock> blockri = real_block_[0]->copy();
    *blockri += *imag_block_[0];

    tmpi = blockri->form_vec(dri);
    *tmpi -= *tmpr;
    *tmpi -= *tmp0;
    *tmpr -= *tmp0;
  }

  auto outr = make_shared<VectorB>(cnaux_);
  auto outi = make_shared<VectorB>(cnaux_);
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
