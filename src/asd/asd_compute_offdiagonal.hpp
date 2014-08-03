//
// BAGEL - Parallel electron correlation program.
// Filename: asd_compute_offdiagonal.hpp
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifdef ASD_HEADERS

#ifndef BAGEL_ASD_COMPUTE_OFFDIAGONAL_H
#define BAGEL_ASD_COMPUTE_OFFDIAGONAL_H

namespace {
  template<typename T>
  void transpose_call(std::shared_ptr<T>& o) { assert(false); }
  template<>
  void transpose_call(std::shared_ptr<Matrix>& o) { o = o->transpose(); }
  template<>
  void transpose_call(std::shared_ptr<RDM<2>>& o) { /* doing nothing */ }
}

template <bool _N, typename return_type>
std::shared_ptr<return_type> ASD_base::couple_blocks(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) const {

  Coupling term_type = coupling_type(AB, ApBp);

  const DimerSubspace_base* space1 = &AB;
  const DimerSubspace_base* space2 = &ApBp;

  bool flip = (static_cast<int>(term_type) < 0);

  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    std::swap(space1,space2);
  }

  std::shared_ptr<return_type> out;

  std::array<MonomerKey,4> keys {{space1->template monomerkey<0>(), space1->template monomerkey<1>(), space2->template monomerkey<0>(), space2->template monomerkey<1>()}};

  switch(term_type) {
    case Coupling::none :
      out = nullptr; break;
    case Coupling::diagonal :
      out = compute_inter_2e<_N>(keys); break;
    case Coupling::aET :
      out = compute_aET<_N>(keys); break;
    case Coupling::bET :
      out = compute_bET<_N>(keys); break;
    case Coupling::abFlip :
      out = compute_abFlip<_N>(keys); break;
    case Coupling::abET :
      out = compute_abET<_N>(keys); break;
    case Coupling::aaET :
      out = compute_aaET<_N>(keys); break;
    case Coupling::bbET :
      out = compute_bbET<_N>(keys); break;
    default :
      throw std::logic_error("Asking for a coupling type that has not been written.");
  }

  /* if we are computing the Hamiltonian and flip = true, then we tranpose the output (see above) */
  if (flip) transpose_call(out);

  return out;
}

#endif

#endif
