//
// BAGEL - Parallel electron correlation program.
// Filename: dfhalfcomplex.cc
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


#include <src/rel/dfhalfcomplex.h>

using namespace std;
using namespace bagel;

DFHalfComplex::DFHalfComplex(const shared_ptr<const DFDist> df, shared_ptr<const Matrix> rcoeff, shared_ptr<const Matrix> icoeff, 
                             const bool swap, pair<const int, const int> coord) : coord_(coord) {

  shared_ptr<DFHalfDist> rhalfbj;
  shared_ptr<DFHalfDist> ihalfbj;

  if (swap == true) {
    rhalfbj = df->compute_half_transform_swap(rcoeff);
    ihalfbj = df->compute_half_transform_swap(icoeff); 
  } else {
    rhalfbj = df->compute_half_transform(rcoeff);
    ihalfbj = df->compute_half_transform(icoeff); 
  }

  dfdata_[0] = rhalfbj->apply_J();
  dfdata_[1] = ihalfbj->apply_J();
}

#if 0
pair<pair<const int, const int>, shared_ptr<ZMatrix> > DFHalfComplex::form_2index_complex(shared_ptr<const DFHalfComplex> dfc) {
  data_[0]->zero();
  data_[1]->zero();
  *data_[0] += *dfdata_[0]->form_2index(dfc->get_real(), 1.0);
  *data_[0] += *dfdata_[1]->form_2index(dfc->get_imag(), -1.0);
  *data_[1] += *dfdata_[0]->form_2index(dfc->get_imag(), 1.0);
  *data_[1] += *dfdata_[1]->form_2index(dfc->get_real(), 1.0);

}
#endif
