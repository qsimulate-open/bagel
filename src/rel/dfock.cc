//
// BAGEL - Parallel electron correlation program.
// Filename: dfock.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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
#include <src/rel/dfock.h>
#include <src/util/matrix.h>
#include <src/rel/smalleribatch.h>

using namespace std;
using namespace bagel;

void DFock::two_electron_part(const std::shared_ptr<const ZMatrix> ocoeff, const bool rhf, const double scale_exchange) {

  if (!rhf) throw logic_error("DFock::two_electron_part() is not implemented for non RHF cases");

//large part
  complex<double> imag (0.0,1.0);
  std::shared_ptr<const DFDist> df = geom_->df();
//std::shared_ptr<const DFDist> dfs = geom_->form_fit<DFDist_ints<SmallERIBatch> >(1.0e-8, false); // TODO thresh should be controlled from the input deck

#if 0
 dfs[0] = shared_ptr<DFDist_ints>(new DFDist_ints(geom_->nbasis(), geom_->naux(), batch[0]))
 dfs[1] = shared_ptr<DFDist_ints>(new DFDist_ints(batch[1]))
 dfs[2] = shared_ptr<DFDist_ints>(new DFDist_ints(batch[2]))
 dfs[3] = shared_ptr<DFDist_ints>(new DFDist_ints(batch[3]))
#endif

  std::shared_ptr<const Matrix> rocoeff = ocoeff->get_real_part(); 
  std::shared_ptr<const Matrix> iocoeff = ocoeff->get_imag_part(); 

  std::shared_ptr<DFHalfDist> rhalfbj = df->compute_half_transform(rocoeff);
  std::shared_ptr<DFHalfDist> ihalfbj = df->compute_half_transform(iocoeff);
  std::shared_ptr<DFHalfDist> rhalf = rhalfbj->apply_J();
  std::shared_ptr<DFHalfDist> ihalf = ihalfbj->apply_J();

#if 0
  *this += *rhalf->form_2index(rhalf, -1.0*scale_exchange);
  *this += *ihalf->form_2index(ihalf, -1.0*scale_exchange);
  shared_ptr<Matrix> cross = rhalf->form_2index(ihalf, -2.0*scale_exchange);
  // TODO maybe not needed
  cross->symmetrize();
  *this += *cross * imag;

  std::shared_ptr<const Matrix> trocoeff(new Matrix(*rocoeff->transpose()*2.0));
  std::shared_ptr<const Matrix> tiocoeff(new Matrix(*iocoeff->transpose()*2.0));
  *this += *df->compute_Jop(rhalfbj, trocoeff);
  *this -= *df->compute_Jop(ihalfbj, tiocoeff);
  const Matrix = *df->compute_Jop(ihalfbj, trocoeff) + *df->compute_Jop(rhalfbj, tiocoeff);
  ZMatrix xx
  *this += *df->compute_Jop(ihalfbj, trocoeff) * imag;
  *this += *df->compute_Jop(rhalfbj, tiocoeff) * imag;
#endif

}
