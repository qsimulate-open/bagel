//
// Newint - Parallel electron correlation program.
// Filename: dipole.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/prop/dipole.h>
#include <src/osint/dipolebatch.h>
#include <iomanip>
#include <array>

using namespace std;

Dipole::Dipole(shared_ptr<const Geometry> g, shared_ptr<const Matrix1e> d) : geom_(g), den_(d) {

}


Dipole::~Dipole() {

}


array<double,3> Dipole::compute() const {
  array<double,3> out{{0.0, 0.0, 0.0}};
  array<double,3> center = geom_->charge_center();

  const vector<vector<int> > offsets = geom_->offsets();

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
          DipoleBatch dipole(input, center);
          dipole.compute();

          const int dimb1 = input[0]->nbasis(); 
          const int dimb0 = input[1]->nbasis(); 
          const double* dat0 = dipole.data();
          const double* dat1 = dipole.data() + dipole.size_block();
          const double* dat2 = dipole.data() + dipole.size_block()*2;
          for (int i = *offset0; i != dimb0 + *offset0; ++i) {
            for (int j = *offset1; j != dimb1 + *offset1; ++j, ++dat0, ++dat1, ++dat2) {
              out[0] += *dat0 * den_->element(j,i);
              out[1] += *dat1 * den_->element(j,i);
              out[2] += *dat2 * den_->element(j,i);
            }
          }

        }
      }
    }
  }

  cout << "    * Permanent dipole moment: (" << setw(12) << setprecision(6) << out[0] << ", "
                                             << setw(12) << out[1] << ", " << setw(12) << out[2] << ") a.u." << endl;

  return out;
}
