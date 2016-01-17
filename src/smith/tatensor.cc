//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: tatensor.cc
// Copyright (C) 2015 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/tatensor.h>

template class bagel::SMITH::TATensor<double,1>;
template class bagel::SMITH::TATensor<double,2>;
template class bagel::SMITH::TATensor<double,3>;
template class bagel::SMITH::TATensor<double,4>;
template class bagel::SMITH::TATensor<double,5>;
template class bagel::SMITH::TATensor<double,6>;
template class bagel::SMITH::TATensor<double,7>;
template class bagel::SMITH::TATensor<std::complex<double>,1>;
template class bagel::SMITH::TATensor<std::complex<double>,2>;
template class bagel::SMITH::TATensor<std::complex<double>,3>;
template class bagel::SMITH::TATensor<std::complex<double>,4>;
template class bagel::SMITH::TATensor<std::complex<double>,5>;
template class bagel::SMITH::TATensor<std::complex<double>,6>;
template class bagel::SMITH::TATensor<std::complex<double>,7>;

#endif
