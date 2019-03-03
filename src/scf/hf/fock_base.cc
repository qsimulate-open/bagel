//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fock_base.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/scf/hf/fock.h>

using namespace std;
using namespace bagel;

template <typename MatType, class Enable>
Fock_base_<MatType, Enable>::Fock_base_(shared_ptr<const Geometry> geom, shared_ptr<const MatType> previous, std::shared_ptr<const MatType> den, const vector<double>& schwarz)
 : Matrix1e_<MatType>(geom), geom_(geom), previous_(previous), density_(den), schwarz_(schwarz) {

  schwarz_thresh_ = geom->schwarz_thresh();

  init(geom); // zero here
}


template <typename MatType, class Enable>
void Fock_base_<MatType, Enable>::fock_one_electron_part() {

  assert(ndim() == mdim());

  *this += *previous_;

  fill_upper_conjg();
}


template <typename MatType, class Enable>
void Fock_base_<MatType, Enable>::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule>) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb0 = input[1]->nbasis();
  const int dimb1 = input[0]->nbasis();

  for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
    for (int j = offsetb1; j != dimb1 + offsetb1; ++j) {
      element(j, i) = 0.0;
    }
  }
}


template class bagel::Fock_base_<Matrix>;
template class bagel::Fock_base_<ZMatrix>;

BOOST_CLASS_EXPORT_IMPLEMENT(Fock_base)
BOOST_CLASS_EXPORT_IMPLEMENT(Fock_base_London)

