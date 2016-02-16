//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: angmom_london.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#include <src/mat1e/giao/angmom_london.h>
#include <src/integral/compos/complexangmombatch.h>
#include <src/integral/comprys/test_codes/bagel_interface.h>

using namespace std;
using namespace bagel;

AngMom_London::AngMom_London(shared_ptr<const Geometry> g, array<double,3> mc) : geom_(g), mcoord_(mc) {
  assert(geom_->magnetism());
}


array<shared_ptr<ZMatrix>, 3> AngMom_London::compute() const {

  const int nbasis = geom_->nbasis();
  auto mat0 = make_shared<ZMatrix>(nbasis, nbasis);
  auto mat1 = make_shared<ZMatrix>(nbasis, nbasis);
  auto mat2 = make_shared<ZMatrix>(nbasis, nbasis);

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
          ComplexAngMomBatch mom(input, geom_->magnetic_field(), mcoord_);
          mom.compute();
          /**************************/
#if 0
          bool correct = true;
          cout << setprecision(12) << fixed;
          const std::complex<double>* angmomdata = mom.data();
          std::vector<std::pair<std::vector<int>,std::complex<double>>> reference = test::get_comparison_orb_angular (input, geom_->magnetic_field(), mcoord_);
          const size_t bagelsize = mom.size_block();
          const size_t refsize = reference.size() / 3;
          assert(refsize * 3 == reference.size());
          const array<char,3> dim = {{'x', 'y', 'z'}};
          cout << endl;
            cout << endl;
            for (int i = 0; i != refsize; ++i) {
          for (int j = 0; j != 3; ++j) {
              cout << "  " << setw(4) << i << " " << dim[j] << " Bagel " << setw(34) << angmomdata[i+j*bagelsize] << "  Test: ";
              for (int k=0; k!=4; ++k) cout << reference[i+j*refsize].first[k] << " ";
              cout << setw(34) << reference[i+j*refsize].second;
              const complex<double> difference = reference[i+j*refsize].second - angmomdata[i+j*bagelsize];
              cout << "  ...  difference = " << setw(34) << difference << endl;
              if (std::abs(difference) > 1.0e-12) correct = false;
            }
          }
          assert(correct);
#endif
          /*************************/
          const complex<double>* dat0 = mom.data();
          const complex<double>* dat1 = mom.data() + mom.size_block();
          const complex<double>* dat2 = mom.data() + mom.size_block()*2;
          for (int i = *offset0; i != *offset0 + (*b0)->nbasis(); ++i) {
            for (int j = *offset1; j != *offset1 + (*b1)->nbasis(); ++j, ++dat0, ++dat1, ++dat2) {
              mat0->element(j,i) = *dat0;
              mat1->element(j,i) = *dat1;
              mat2->element(j,i) = *dat2;
            }
          }

        }
      }
    }
  }

  return array<shared_ptr<ZMatrix>,3>{{mat0, mat1, mat2}};
}
