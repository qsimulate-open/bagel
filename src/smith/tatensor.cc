//
// BAGEL - Parallel electron correlation program.
// Filename: tatensor.cc
// Copyright (C) 2015 Toru Shiozaki
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
