//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dvec.cc
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

#include <stdexcept>
#include <src/ci/fci/dvec.h>

template <typename DataType>
std::shared_ptr<bagel::Dvector<DataType>> bagel::Dvector<DataType>::extract_state(const std::vector<int> input) const {
  auto out = std::make_shared<bagel::Dvector<DataType>>(det(), input.size());

  for (int i = 0; i != input.size(); ++i)
    std::copy_n(data(input[i])->data(), lenb_*lena_, out->data(i)->data());
    //out->data(i) = data(input[i])->copy();
  return out;
}

template class bagel::Dvector<double>;
template class bagel::Dvector<std::complex<double>>;

