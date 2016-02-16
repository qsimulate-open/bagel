//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: bvrr_optim.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/integral/rys/breitbatch.h>

using namespace std;
using namespace bagel;


// hand written code
void BreitBatch::perform_VRR1() {
  // VRR1 means that angular numbers are both zero
  double t[3]; // tilde
  double s[3]; // tilde,tilde

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const int ii3 = ii*3;
    const double oxpq2 = 0.5 * (xp_[ii]+xq_[ii]) / (xp_[ii]*xq_[ii]);
    const double rho2 = 1.0 / oxpq2 *coeff_[ii]*weights_[ii];
    const double root = 1.0 - roots_[ii];

    for (int i = 0; i != 3; ++i) {
      const double PQ = P_[ii3+i] - Q_[ii3+i];
      t[i] = PQ;
      s[i] = (root*PQ*PQ + oxpq2)*rho2;
    }

    const double fac = root*rho2;

    data_[ii+0*size_block_] = s[0];
    data_[ii+1*size_block_] = t[0]*t[1]*fac;
    data_[ii+2*size_block_] = t[0]*t[2]*fac;
    data_[ii+3*size_block_] = s[1];
    data_[ii+4*size_block_] = t[1]*t[2]*fac;
    data_[ii+5*size_block_] = s[2];

  }
}


void Spin2Batch::perform_VRR1() {
  perform_VRR();
}


