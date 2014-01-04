//
// BAGEL - Parallel electron correlation program.
// Filename: r2batch.cc
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


#include <src/integral/rys/r2batch.h>
#include <src/util/constants.h>
#include <src/integral/rys/r2rootlist.h>

using namespace std;
using namespace bagel;

static constexpr double T_thresh__ = 1.0e-8;

void R2Batch::compute_ssss(const double integral_thresh) {
  screening_size_ = 0;

  const vector<double> exp0 = basisinfo_[0]->exponents();
  const vector<double> exp1 = basisinfo_[1]->exponents();

  int index = 0;
  vector<shared_ptr<const Atom>> atoms = mol_->atoms();

  const double onepi2 = 1.0 / (pi__ * pi__);
  const double sqrtpi = sqrt(pi__);
  for (auto expi0 = exp0.begin(); expi0 != exp0.end(); ++expi0) {
    for (auto expi1 = exp1.begin(); expi1 != exp1.end(); ++expi1) {
      for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++index) {
        double zeta = (*aiter)->ecp(0);
        double Z = (*aiter)->atom_charge();
        const double cxp = *expi0 + *expi1;
        const double socxp = cxp + zeta;
        xp_[index] = cxp;
        const double ab = *expi0 * *expi1;
        const double cxp_inv = 1.0 / cxp;
        const double socxp_inv = 1.0 / socxp;
        P_[index * 3    ] = (basisinfo_[0]->position(0) * *expi0 + basisinfo_[1]->position(0) * *expi1) * cxp_inv;
        P_[index * 3 + 1] = (basisinfo_[0]->position(1) * *expi0 + basisinfo_[1]->position(1) * *expi1) * cxp_inv;
        P_[index * 3 + 2] = (basisinfo_[0]->position(2) * *expi0 + basisinfo_[1]->position(2) * *expi1) * cxp_inv;
        const double Eab = exp(-(AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2]) * (ab * cxp_inv) );
        const double PCx = P_[index * 3    ] - (*aiter)->position(0);
        const double PCy = P_[index * 3 + 1] - (*aiter)->position(1);
        const double PCz = P_[index * 3 + 2] - (*aiter)->position(2);
        const double T = cxp * cxp * socxp_inv * (PCx * PCx + PCy * PCy + PCz * PCz);
        coeff_[index] = 2 * pi__ * sqrtpi * sqrt(socxp_inv) * exp(-cxp * (PCx * PCx + PCy * PCy + PCz * PCz) + T) * Eab;
        const double sqrtt = sqrt(T);
        const double ss = coeff_[index] * pow(4.0 * ab * onepi2, 0.75) * (T > 1.0e-15 ? exp(sqrtt) * inline_dawson(sqrtt) / sqrtt : 1.0);
        if (ss > integral_thresh) {
          T_[index] = T;
          screening_[screening_size_] = index;
          ++screening_size_;
        } else {
          T_[index] = -1.0;
          coeff_[index] = 0.0;
        }
      }
    }
  }
}

void R2Batch::root_weight(const int ps) {
  if (amax_ + cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int i = screening_[j];
      if (T_[i] < T_thresh__) {
        weights_[i] = 1.0;
      } else {
        const double sqrtt = sqrt(T_[i]);
        const double dawsonsqt = inline_dawson(sqrtt);
        weights_[i] = dawsonsqt / sqrtt;
      }
    }
  } else {
    r2root__.root(rank_, T_, roots_, weights_, ps);
  }
}

