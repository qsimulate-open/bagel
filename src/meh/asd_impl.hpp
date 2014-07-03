//
// BAGEL - Parallel electron correlation program.
// Filename: asd_impl.hpp
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifdef MEH_HEADERS

#ifndef BAGEL_MEH_ASD_IMPL_H
#define BAGEL_MEH_ASD_IMPL_H

template <>
template <class VecType>
std::shared_ptr<Matrix> asd::ASD_impl<true>::compute_diagonal_block(MultiExcitonHamiltonian<VecType>* me, DSb<VecType>& subspace) {
  const double core = me->dimer_->sref()->geom()->nuclear_repulsion() + me->jop_->core_energy();

  auto out = me->compute_intra(subspace, me->jop_, core);
  std::array<MonomerKey,4> keys {{ subspace.template monomerkey<0>(), subspace.template monomerkey<1>(),
                                   subspace.template monomerkey<0>(), subspace.template monomerkey<1>() }};
  *out += *me->template compute_inter_2e<true>(keys);

  return out;
}

#endif
#endif
