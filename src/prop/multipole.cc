//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multipole.cc
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


#include <src/prop/multipole.h>
#include <src/integral/os/mmbatch.h>

using namespace std;
using namespace bagel;

Multipole::Multipole(shared_ptr<const Geometry> g, shared_ptr<const Matrix> d, const int rank, const string jobn)
  : geom_(g), den_(d), rank_(rank), jobname_(jobn) {

}


vector<double> Multipole::compute() const {
  vector<double> out;

  // TODO perhaps we could reduce operation by a factor of 2
  auto o0 = geom_->offsets().begin();
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++o0) {
    auto o1 = geom_->offsets().begin();
    for (auto a1 = geom_->atoms().begin(); a1 != geom_->atoms().end(); ++a1, ++o1) {

      auto offset0 = o0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++offset0) {
        auto offset1 = o1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++offset1) {

          array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          MMBatch mpole(input, geom_, rank_);
          mpole.compute();

          if (out.empty())
            out.resize(mpole.num_blocks());

          const int dimb1 = input[0]->nbasis();
          const int dimb0 = input[1]->nbasis();
          vector<const double*> dat(mpole.num_blocks());
          for (int i = 0; i != mpole.num_blocks(); ++i)
            dat[i] = mpole.data() + mpole.size_block()*i;

          for (int i = *offset0; i != dimb0 + *offset0; ++i) {
            for (int j = *offset1; j != dimb1 + *offset1; ++j) {
              for (int k = 0; k != mpole.num_blocks(); ++k) {
                out[k] += *dat[k]++ * den_->element(j,i);
              }
            }
          }

        }
      }
    }
  }

  cout << "    * Permanent dipole moment:" << (jobname_.empty() ? "" : " " + jobname_) << endl;
  cout << "           (" << setw(12) << setprecision(6) << out[0] << ", " << setw(12) << out[1] << ", " << setw(12) << out[2] << ") a.u." << endl << endl;

  if (rank_ == 2) {
    // we need to add nuclear contribution
    // \hat xx = -\sum_i x_i^2 + \sum_N Z_N X_N^2 (same as Molpro's convention).
    array<double,6> sm = geom_->quadrupole();
    for (int i = 0; i != 6; ++i)
      out[i+3] = sm[i] - out[i+3];
    cout << "    * Permanent second moment (around the center of charge):" << (jobname_.empty() ? "" : " " + jobname_) << endl;
    cout << "          Qxx " << setw(12) << setprecision(6) << out[3] << endl;
    cout << "          Qxy " << setw(12) << setprecision(6) << out[4] << endl;
    cout << "          Qyy " << setw(12) << setprecision(6) << out[5] << endl;
    cout << "          Qxz " << setw(12) << setprecision(6) << out[6] << endl;
    cout << "          Qyz " << setw(12) << setprecision(6) << out[7] << endl;
    cout << "          Qzz " << setw(12) << setprecision(6) << out[8] << endl << endl;

  } else if (rank_ >= 3) {
    throw logic_error("higher-order multipole integrals are implemented, but post-processing is missing");
  }

  return out;
}
