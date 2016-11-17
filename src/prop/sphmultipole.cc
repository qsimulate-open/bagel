//
// BAGEL - Parallel electron correlation program.
// Filename: sphmultipole.cc
// Copyright (C) 2016 Toru Shiozaki
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


#include <src/prop/sphmultipole.h>

using namespace std;
using namespace bagel;

SphMultipole::SphMultipole(shared_ptr<const Geometry> g, shared_ptr<const Matrix> d, const int rank)
  : geom_(g), density_(d), rank_(rank) { }

vector<complex<double>> SphMultipole::compute() const { // slow
  
  const size_t size = (rank_+1)*(rank_+1);
  vector<complex<double>> out(size);

  auto o0 = geom_->offsets().begin();
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++o0) {
    auto o1 = geom_->offsets().begin();
    for (auto a1 = geom_->atoms().begin(); a1 != geom_->atoms().end(); ++a1, ++o1) {

      auto offset0 = o0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++offset0) {
        auto offset1 = o1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++offset1) {

          array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          MultipoleBatch mpole(input, geom_->charge_center(), rank_);
          mpole.compute();
          assert(mpole.num_blocks() == size);

          const int dimb1 = input[0]->nbasis();
          const int dimb0 = input[1]->nbasis();
          vector<const complex<double>*> dat(mpole.num_blocks());
          for (int i = 0; i != mpole.num_blocks(); ++i)
            dat[i] = mpole.data() + mpole.size_block()*i;

          for (int i = *offset0; i != dimb0 + *offset0; ++i) {
            for (int j = *offset1; j != dimb1 + *offset1; ++j) {
              for (int k = 0; k != mpole.num_blocks(); ++k) {
                out[k] += *dat[k]++ * density_->element(j,i);
              }
            }
          }

        }
      }
    }
  }

  cout << "    * Centre of charge:      (" << geom_->charge_center()[0] << ", " << setw(12)
                                           << geom_->charge_center()[1] << ", " << setw(12)
                                           << geom_->charge_center()[2] << ")" << endl;
  cout << "    * Permanent dipole moment:" << endl;
  cout << "           (" << setw(12) << setprecision(6) << out[1] << ", " << setw(12) << out[3] << ", " << setw(12) << out[2] << ") a.u." << endl << endl;

  if (rank_ >= 2) {
    cout << "    * Permanent quadrupole moment:" << endl;
    cout << "          Q2-2 " << setw(11) << setprecision(6) << out[4] << endl;
    cout << "          Q2-1 " << setw(11) << setprecision(6) << out[5] << endl;
    cout << "          Q20  " << setw(11) << setprecision(6) << out[6] << endl;
    cout << "          Q21  " << setw(11) << setprecision(6) << out[7] << endl;
    cout << "          Q22  " << setw(11) << setprecision(6) << out[8] << endl;
   }

  return out;
}
