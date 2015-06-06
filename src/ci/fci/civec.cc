//
// BAGEL - Parallel electron correlation program.
// Filename: civec.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/ci/fci/civec.h>
#include <src/util/math/algo.h>

template class bagel::Civector<double>;
template class bagel::Civector<std::complex<double>>;

using namespace std;
using namespace bagel;

template<>
void Civector<double>::spin_decontaminate(const double thresh) {
  const int nspin = det_->nspin();
  const int max_spin = det_->nelea() + det_->neleb();

  const double expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<Civec> S2 = spin();

  int k = nspin + 2;
  while( fabs(dot_product(*S2) - expectation) > thresh ) {
    if ( k > max_spin ) throw runtime_error("Spin decontamination failed.");

    const double factor = -4.0/(static_cast<double>(k*(k+2)));
    ax_plus_y(factor, *S2);

    const double norm = this->norm();
    const double rescale = (norm*norm > 1.0e-60) ? 1.0/norm : 0.0;
    scale(rescale);

    S2 = spin();

    k += 2;
  }
}


template<>
void CASDvec::apply_and_fill(shared_ptr<const CASDvec> source_dvec, const int orbital, const bool action, const bool spin) {
  shared_ptr<const DetType> source_det = source_dvec->det();
  shared_ptr<const DetType> target_det = this->det();

  this->zero();

  const int source_lenb = source_det->lenb();
  const int target_lenb = target_det->lenb();
  const int target_lena = target_det->lena();

  if (spin) {
    for (size_t ivec = 0; ivec != this->ij(); ++ivec) {
      double* target_base = this->data(ivec)->data();
      const double* source_base = source_dvec->data(ivec)->data();
      for (auto& iter : (action ? source_det->phiupa(orbital) : source_det->phidowna(orbital))) {
        const double sign = static_cast<double>(iter.sign);
        double* target = target_base + target_lenb * iter.target;
        const double* source = source_base + source_lenb * iter.source;
        blas::ax_plus_y_n(sign, source, target_lenb, target);
      }
    }
  } else {
    for (size_t ivec = 0; ivec != this->ij(); ++ivec) {
      for (int i = 0; i < target_lena; ++i) {
        double* target_base = this->data(ivec)->element_ptr(0,i);
        const double* source_base = source_dvec->data(ivec)->element_ptr(0,i);
        for (auto& iter : (action ? source_det->phiupb(orbital) : source_det->phidownb(orbital))) {
          const double sign = static_cast<double>(iter.sign);
          target_base[iter.target] += sign * source_base[iter.source];
        }
      }
    }
  }
}

